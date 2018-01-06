#line 1 "dlamtsqr.f"
/* dlamtsqr.f -- translated by f2c (version 20100827).
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

#line 1 "dlamtsqr.f"
/* Table of constant values */

static integer c__0 = 0;


/*  Definition: */
/*  =========== */

/*      SUBROUTINE DLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, */
/*     $                     LDT, C, LDC, WORK, LWORK, INFO ) */


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
/* >      DLAMTSQR overwrites the general real M-by-N matrix C with */
/* > */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* >      where Q is a real orthogonal matrix defined as the product */
/* >      of blocked elementary reflectors computed by tall skinny */
/* >      QR factorization (DLATSQR) */
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
/* >          The number of columns of the matrix C. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          N >= K >= 0; */
/* > */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* >          MB is INTEGER */
/* >          The block size to be used in the blocked QR. */
/* >          MB > N. (must be the same as DLATSQR) */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The column block size to be used in the blocked QR. */
/* >          N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          blockedelementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by DLATSQR in the first k columns of */
/* >          its array argument A. */
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
/* >          ( N * Number of blocks(CEIL(M-K/MB-K)), */
/* >          The blocked upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks.  See below */
/* >          for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
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
/* > */
/* > \endverbatim */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* > */
/* >          If SIDE = 'L', LWORK >= max(1,N)*NB; */
/* >          if SIDE = 'R', LWORK >= max(1,MB)*NB. */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > */
/* > \endverbatim */
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
/* > Tall-Skinny QR (TSQR) performs QR by a sequence of orthogonal transformations, */
/* > representing Q as a product of other orthogonal matrices */
/* >   Q = Q(1) * Q(2) * . . . * Q(k) */
/* > where each Q(i) zeros out subdiagonal entries of a block of MB rows of A: */
/* >   Q(1) zeros out the subdiagonal entries of rows 1:MB of A */
/* >   Q(2) zeros out the bottom MB-N rows of rows [1:N,MB+1:2*MB-N] of A */
/* >   Q(3) zeros out the bottom MB-N rows of rows [1:N,2*MB-N+1:3*MB-2*N] of A */
/* >   . . . */
/* > */
/* > Q(1) is computed by GEQRT, which represents Q(1) by Householder vectors */
/* > stored under the diagonal of rows 1:MB of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,1:N). */
/* > For more information see Further Details in GEQRT. */
/* > */
/* > Q(i) for i>1 is computed by TPQRT, which represents Q(i) by Householder vectors */
/* > stored in rows [(i-1)*(MB-N)+N+1:i*(MB-N)+N] of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,(i-1)*N+1:i*N). */
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
/* Subroutine */ int dlamtsqr_(char *side, char *trans, integer *m, integer *
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
    extern /* Subroutine */ int dgemqrt_(char *, char *, integer *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), dtpmqrt_(char *, char *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);


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

#line 229 "dlamtsqr.f"
    /* Parameter adjustments */
#line 229 "dlamtsqr.f"
    a_dim1 = *lda;
#line 229 "dlamtsqr.f"
    a_offset = 1 + a_dim1;
#line 229 "dlamtsqr.f"
    a -= a_offset;
#line 229 "dlamtsqr.f"
    t_dim1 = *ldt;
#line 229 "dlamtsqr.f"
    t_offset = 1 + t_dim1;
#line 229 "dlamtsqr.f"
    t -= t_offset;
#line 229 "dlamtsqr.f"
    c_dim1 = *ldc;
#line 229 "dlamtsqr.f"
    c_offset = 1 + c_dim1;
#line 229 "dlamtsqr.f"
    c__ -= c_offset;
#line 229 "dlamtsqr.f"
    --work;
#line 229 "dlamtsqr.f"

#line 229 "dlamtsqr.f"
    /* Function Body */
#line 229 "dlamtsqr.f"
    lquery = *lwork < 0;
#line 230 "dlamtsqr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 231 "dlamtsqr.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 232 "dlamtsqr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 233 "dlamtsqr.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 234 "dlamtsqr.f"
    if (left) {
#line 235 "dlamtsqr.f"
	lw = *n * *nb;
#line 236 "dlamtsqr.f"
    } else {
#line 237 "dlamtsqr.f"
	lw = *mb * *nb;
#line 238 "dlamtsqr.f"
    }

#line 240 "dlamtsqr.f"
    *info = 0;
#line 241 "dlamtsqr.f"
    if (! left && ! right) {
#line 242 "dlamtsqr.f"
	*info = -1;
#line 243 "dlamtsqr.f"
    } else if (! tran && ! notran) {
#line 244 "dlamtsqr.f"
	*info = -2;
#line 245 "dlamtsqr.f"
    } else if (*m < 0) {
#line 246 "dlamtsqr.f"
	*info = -3;
#line 247 "dlamtsqr.f"
    } else if (*n < 0) {
#line 248 "dlamtsqr.f"
	*info = -4;
#line 249 "dlamtsqr.f"
    } else if (*k < 0) {
#line 250 "dlamtsqr.f"
	*info = -5;
#line 251 "dlamtsqr.f"
    } else if (*lda < max(1,*k)) {
#line 252 "dlamtsqr.f"
	*info = -9;
#line 253 "dlamtsqr.f"
    } else if (*ldt < max(1,*nb)) {
#line 254 "dlamtsqr.f"
	*info = -11;
#line 255 "dlamtsqr.f"
    } else if (*ldc < max(1,*m)) {
#line 256 "dlamtsqr.f"
	*info = -13;
#line 257 "dlamtsqr.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 258 "dlamtsqr.f"
	*info = -15;
#line 259 "dlamtsqr.f"
    }

/*     Determine the block size if it is tall skinny or short and wide */

#line 263 "dlamtsqr.f"
    if (*info == 0) {
#line 264 "dlamtsqr.f"
	work[1] = (doublereal) lw;
#line 265 "dlamtsqr.f"
    }

#line 267 "dlamtsqr.f"
    if (*info != 0) {
#line 268 "dlamtsqr.f"
	i__1 = -(*info);
#line 268 "dlamtsqr.f"
	xerbla_("DLAMTSQR", &i__1, (ftnlen)8);
#line 269 "dlamtsqr.f"
	return 0;
#line 270 "dlamtsqr.f"
    } else if (lquery) {
#line 271 "dlamtsqr.f"
	return 0;
#line 272 "dlamtsqr.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 276 "dlamtsqr.f"
    i__1 = min(*m,*n);
#line 276 "dlamtsqr.f"
    if (min(i__1,*k) == 0) {
#line 277 "dlamtsqr.f"
	return 0;
#line 278 "dlamtsqr.f"
    }

/* Computing MAX */
#line 280 "dlamtsqr.f"
    i__1 = max(*m,*n);
#line 280 "dlamtsqr.f"
    if (*mb <= *k || *mb >= max(i__1,*k)) {
#line 281 "dlamtsqr.f"
	dgemqrt_(side, trans, m, n, k, nb, &a[a_offset], lda, &t[t_offset], 
		ldt, &c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)
		1);
#line 283 "dlamtsqr.f"
	return 0;
#line 284 "dlamtsqr.f"
    }

#line 286 "dlamtsqr.f"
    if (left && notran) {

/*         Multiply Q to the last block of C */

#line 290 "dlamtsqr.f"
	kk = (*m - *k) % (*mb - *k);
#line 291 "dlamtsqr.f"
	ctr = (*m - *k) / (*mb - *k);
#line 292 "dlamtsqr.f"
	if (kk > 0) {
#line 293 "dlamtsqr.f"
	    ii = *m - kk + 1;
#line 294 "dlamtsqr.f"
	    dtpmqrt_("L", "N", &kk, n, k, &c__0, nb, &a[ii + a_dim1], lda, &t[
		    (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, 
		    &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 297 "dlamtsqr.f"
	} else {
#line 298 "dlamtsqr.f"
	    ii = *m + 1;
#line 299 "dlamtsqr.f"
	}

#line 301 "dlamtsqr.f"
	i__1 = *mb + 1;
#line 301 "dlamtsqr.f"
	i__2 = -(*mb - *k);
#line 301 "dlamtsqr.f"
	for (i__ = ii - (*mb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (I:I+MB,1:N) */

#line 305 "dlamtsqr.f"
	    --ctr;
#line 306 "dlamtsqr.f"
	    i__3 = *mb - *k;
#line 306 "dlamtsqr.f"
	    dtpmqrt_("L", "N", &i__3, n, k, &c__0, nb, &a[i__ + a_dim1], lda, 
		    &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);

#line 310 "dlamtsqr.f"
	}

/*         Multiply Q to the first block of C (1:MB,1:N) */

#line 314 "dlamtsqr.f"
	dgemqrt_("L", "N", mb, n, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 317 "dlamtsqr.f"
    } else if (left && tran) {

/*         Multiply Q to the first block of C */

#line 321 "dlamtsqr.f"
	kk = (*m - *k) % (*mb - *k);
#line 322 "dlamtsqr.f"
	ii = *m - kk + 1;
#line 323 "dlamtsqr.f"
	ctr = 1;
#line 324 "dlamtsqr.f"
	dgemqrt_("L", "T", mb, n, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 327 "dlamtsqr.f"
	i__2 = ii - *mb + *k;
#line 327 "dlamtsqr.f"
	i__1 = *mb - *k;
#line 327 "dlamtsqr.f"
	for (i__ = *mb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (I:I+MB,1:N) */

#line 331 "dlamtsqr.f"
	    i__3 = *mb - *k;
#line 331 "dlamtsqr.f"
	    dtpmqrt_("L", "T", &i__3, n, k, &c__0, nb, &a[i__ + a_dim1], lda, 
		    &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 334 "dlamtsqr.f"
	    ++ctr;

#line 336 "dlamtsqr.f"
	}
#line 337 "dlamtsqr.f"
	if (ii <= *m) {

/*         Multiply Q to the last block of C */

#line 341 "dlamtsqr.f"
	    dtpmqrt_("L", "T", &kk, n, k, &c__0, nb, &a[ii + a_dim1], lda, &t[
		    (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, 
		    &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);

#line 345 "dlamtsqr.f"
	}

#line 347 "dlamtsqr.f"
    } else if (right && tran) {

/*         Multiply Q to the last block of C */

#line 351 "dlamtsqr.f"
	kk = (*n - *k) % (*mb - *k);
#line 352 "dlamtsqr.f"
	ctr = (*n - *k) / (*mb - *k);
#line 353 "dlamtsqr.f"
	if (kk > 0) {
#line 354 "dlamtsqr.f"
	    ii = *n - kk + 1;
#line 355 "dlamtsqr.f"
	    dtpmqrt_("R", "T", m, &kk, k, &c__0, nb, &a[ii + a_dim1], lda, &t[
		    (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, 
		    &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 358 "dlamtsqr.f"
	} else {
#line 359 "dlamtsqr.f"
	    ii = *n + 1;
#line 360 "dlamtsqr.f"
	}

#line 362 "dlamtsqr.f"
	i__1 = *mb + 1;
#line 362 "dlamtsqr.f"
	i__2 = -(*mb - *k);
#line 362 "dlamtsqr.f"
	for (i__ = ii - (*mb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 366 "dlamtsqr.f"
	    --ctr;
#line 367 "dlamtsqr.f"
	    i__3 = *mb - *k;
#line 367 "dlamtsqr.f"
	    dtpmqrt_("R", "T", m, &i__3, k, &c__0, nb, &a[i__ + a_dim1], lda, 
		    &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);

#line 371 "dlamtsqr.f"
	}

/*         Multiply Q to the first block of C (1:M,1:MB) */

#line 375 "dlamtsqr.f"
	dgemqrt_("R", "T", m, mb, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 378 "dlamtsqr.f"
    } else if (right && notran) {

/*         Multiply Q to the first block of C */

#line 382 "dlamtsqr.f"
	kk = (*n - *k) % (*mb - *k);
#line 383 "dlamtsqr.f"
	ii = *n - kk + 1;
#line 384 "dlamtsqr.f"
	ctr = 1;
#line 385 "dlamtsqr.f"
	dgemqrt_("R", "N", m, mb, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 388 "dlamtsqr.f"
	i__2 = ii - *mb + *k;
#line 388 "dlamtsqr.f"
	i__1 = *mb - *k;
#line 388 "dlamtsqr.f"
	for (i__ = *mb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 392 "dlamtsqr.f"
	    i__3 = *mb - *k;
#line 392 "dlamtsqr.f"
	    dtpmqrt_("R", "N", m, &i__3, k, &c__0, nb, &a[i__ + a_dim1], lda, 
		    &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 395 "dlamtsqr.f"
	    ++ctr;

#line 397 "dlamtsqr.f"
	}
#line 398 "dlamtsqr.f"
	if (ii <= *n) {

/*         Multiply Q to the last block of C */

#line 402 "dlamtsqr.f"
	    dtpmqrt_("R", "N", m, &kk, k, &c__0, nb, &a[ii + a_dim1], lda, &t[
		    (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, 
		    &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);

#line 406 "dlamtsqr.f"
	}

#line 408 "dlamtsqr.f"
    }

#line 410 "dlamtsqr.f"
    work[1] = (doublereal) lw;
#line 411 "dlamtsqr.f"
    return 0;

/*     End of DLAMTSQR */

} /* dlamtsqr_ */

