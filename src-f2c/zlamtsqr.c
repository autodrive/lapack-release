#line 1 "zlamtsqr.f"
/* zlamtsqr.f -- translated by f2c (version 20100827).
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

#line 1 "zlamtsqr.f"
/* Table of constant values */

static integer c__0 = 0;


/*  Definition: */
/*  =========== */

/*      SUBROUTINE ZLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, */
/*     $                     LDT, C, LDC, WORK, LWORK, INFO ) */


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
/* >      ZLAMTSQR overwrites the general complex M-by-N matrix C with */
/* > */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* >      where Q is a real orthogonal matrix defined as the product */
/* >      of blocked elementary reflectors computed by tall skinny */
/* >      QR factorization (ZLATSQR) */
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
/* >          = 'C':  Conjugate Transpose, apply Q**H. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,K) */
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
/* >          T is COMPLEX*16 array, dimension */
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
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
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
/* >         (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* Subroutine */ int zlamtsqr_(char *side, char *trans, integer *m, integer *
	n, integer *k, integer *mb, integer *nb, doublecomplex *a, integer *
	lda, doublecomplex *t, integer *ldt, doublecomplex *c__, integer *ldc,
	 doublecomplex *work, integer *lwork, integer *info, ftnlen side_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    extern /* Subroutine */ int ztpmqrt_(char *, char *, integer *, integer *,
	     integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    static integer i__, ii, kk, lw, ctr;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran, lquery;
    extern /* Subroutine */ int zgemqrt_(char *, char *, integer *, integer *,
	     integer *, integer *, doublecomplex *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, doublecomplex *, integer 
	    *, ftnlen, ftnlen);


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

#line 229 "zlamtsqr.f"
    /* Parameter adjustments */
#line 229 "zlamtsqr.f"
    a_dim1 = *lda;
#line 229 "zlamtsqr.f"
    a_offset = 1 + a_dim1;
#line 229 "zlamtsqr.f"
    a -= a_offset;
#line 229 "zlamtsqr.f"
    t_dim1 = *ldt;
#line 229 "zlamtsqr.f"
    t_offset = 1 + t_dim1;
#line 229 "zlamtsqr.f"
    t -= t_offset;
#line 229 "zlamtsqr.f"
    c_dim1 = *ldc;
#line 229 "zlamtsqr.f"
    c_offset = 1 + c_dim1;
#line 229 "zlamtsqr.f"
    c__ -= c_offset;
#line 229 "zlamtsqr.f"
    --work;
#line 229 "zlamtsqr.f"

#line 229 "zlamtsqr.f"
    /* Function Body */
#line 229 "zlamtsqr.f"
    lquery = *lwork < 0;
#line 230 "zlamtsqr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 231 "zlamtsqr.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 232 "zlamtsqr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 233 "zlamtsqr.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 234 "zlamtsqr.f"
    if (left) {
#line 235 "zlamtsqr.f"
	lw = *n * *nb;
#line 236 "zlamtsqr.f"
    } else {
#line 237 "zlamtsqr.f"
	lw = *m * *nb;
#line 238 "zlamtsqr.f"
    }

#line 240 "zlamtsqr.f"
    *info = 0;
#line 241 "zlamtsqr.f"
    if (! left && ! right) {
#line 242 "zlamtsqr.f"
	*info = -1;
#line 243 "zlamtsqr.f"
    } else if (! tran && ! notran) {
#line 244 "zlamtsqr.f"
	*info = -2;
#line 245 "zlamtsqr.f"
    } else if (*m < 0) {
#line 246 "zlamtsqr.f"
	*info = -3;
#line 247 "zlamtsqr.f"
    } else if (*n < 0) {
#line 248 "zlamtsqr.f"
	*info = -4;
#line 249 "zlamtsqr.f"
    } else if (*k < 0) {
#line 250 "zlamtsqr.f"
	*info = -5;
#line 251 "zlamtsqr.f"
    } else if (*lda < max(1,*k)) {
#line 252 "zlamtsqr.f"
	*info = -9;
#line 253 "zlamtsqr.f"
    } else if (*ldt < max(1,*nb)) {
#line 254 "zlamtsqr.f"
	*info = -11;
#line 255 "zlamtsqr.f"
    } else if (*ldc < max(1,*m)) {
#line 256 "zlamtsqr.f"
	*info = -13;
#line 257 "zlamtsqr.f"
    } else if (*lwork < max(1,lw) && ! lquery) {
#line 258 "zlamtsqr.f"
	*info = -15;
#line 259 "zlamtsqr.f"
    }

/*     Determine the block size if it is tall skinny or short and wide */

#line 263 "zlamtsqr.f"
    if (*info == 0) {
#line 264 "zlamtsqr.f"
	work[1].r = (doublereal) lw, work[1].i = 0.;
#line 265 "zlamtsqr.f"
    }

#line 267 "zlamtsqr.f"
    if (*info != 0) {
#line 268 "zlamtsqr.f"
	i__1 = -(*info);
#line 268 "zlamtsqr.f"
	xerbla_("ZLAMTSQR", &i__1, (ftnlen)8);
#line 269 "zlamtsqr.f"
	return 0;
#line 270 "zlamtsqr.f"
    } else if (lquery) {
#line 271 "zlamtsqr.f"
	return 0;
#line 272 "zlamtsqr.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 276 "zlamtsqr.f"
    i__1 = min(*m,*n);
#line 276 "zlamtsqr.f"
    if (min(i__1,*k) == 0) {
#line 277 "zlamtsqr.f"
	return 0;
#line 278 "zlamtsqr.f"
    }

/* Computing MAX */
#line 280 "zlamtsqr.f"
    i__1 = max(*m,*n);
#line 280 "zlamtsqr.f"
    if (*mb <= *k || *mb >= max(i__1,*k)) {
#line 281 "zlamtsqr.f"
	zgemqrt_(side, trans, m, n, k, nb, &a[a_offset], lda, &t[t_offset], 
		ldt, &c__[c_offset], ldc, &work[1], info, (ftnlen)1, (ftnlen)
		1);
#line 283 "zlamtsqr.f"
	return 0;
#line 284 "zlamtsqr.f"
    }

#line 286 "zlamtsqr.f"
    if (left && notran) {

/*         Multiply Q to the last block of C */

#line 290 "zlamtsqr.f"
	kk = (*m - *k) % (*mb - *k);
#line 291 "zlamtsqr.f"
	ctr = (*m - *k) / (*mb - *k);
#line 292 "zlamtsqr.f"
	if (kk > 0) {
#line 293 "zlamtsqr.f"
	    ii = *m - kk + 1;
#line 294 "zlamtsqr.f"
	    ztpmqrt_("L", "N", &kk, n, k, &c__0, nb, &a[ii + a_dim1], lda, &t[
		    (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, 
		    &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 297 "zlamtsqr.f"
	} else {
#line 298 "zlamtsqr.f"
	    ii = *m + 1;
#line 299 "zlamtsqr.f"
	}

#line 301 "zlamtsqr.f"
	i__1 = *mb + 1;
#line 301 "zlamtsqr.f"
	i__2 = -(*mb - *k);
#line 301 "zlamtsqr.f"
	for (i__ = ii - (*mb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (I:I+MB,1:N) */

#line 305 "zlamtsqr.f"
	    --ctr;
#line 306 "zlamtsqr.f"
	    i__3 = *mb - *k;
#line 306 "zlamtsqr.f"
	    ztpmqrt_("L", "N", &i__3, n, k, &c__0, nb, &a[i__ + a_dim1], lda, 
		    &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 310 "zlamtsqr.f"
	}

/*         Multiply Q to the first block of C (1:MB,1:N) */

#line 314 "zlamtsqr.f"
	zgemqrt_("L", "N", mb, n, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 317 "zlamtsqr.f"
    } else if (left && tran) {

/*         Multiply Q to the first block of C */

#line 321 "zlamtsqr.f"
	kk = (*m - *k) % (*mb - *k);
#line 322 "zlamtsqr.f"
	ii = *m - kk + 1;
#line 323 "zlamtsqr.f"
	ctr = 1;
#line 324 "zlamtsqr.f"
	zgemqrt_("L", "C", mb, n, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 327 "zlamtsqr.f"
	i__2 = ii - *mb + *k;
#line 327 "zlamtsqr.f"
	i__1 = *mb - *k;
#line 327 "zlamtsqr.f"
	for (i__ = *mb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (I:I+MB,1:N) */

#line 331 "zlamtsqr.f"
	    i__3 = *mb - *k;
#line 331 "zlamtsqr.f"
	    ztpmqrt_("L", "C", &i__3, n, k, &c__0, nb, &a[i__ + a_dim1], lda, 
		    &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[i__ + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 334 "zlamtsqr.f"
	    ++ctr;

#line 336 "zlamtsqr.f"
	}
#line 337 "zlamtsqr.f"
	if (ii <= *m) {

/*         Multiply Q to the last block of C */

#line 341 "zlamtsqr.f"
	    ztpmqrt_("L", "C", &kk, n, k, &c__0, nb, &a[ii + a_dim1], lda, &t[
		    (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, 
		    &c__[ii + c_dim1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);

#line 345 "zlamtsqr.f"
	}

#line 347 "zlamtsqr.f"
    } else if (right && tran) {

/*         Multiply Q to the last block of C */

#line 351 "zlamtsqr.f"
	kk = (*n - *k) % (*mb - *k);
#line 352 "zlamtsqr.f"
	ctr = (*n - *k) / (*mb - *k);
#line 353 "zlamtsqr.f"
	if (kk > 0) {
#line 354 "zlamtsqr.f"
	    ii = *n - kk + 1;
#line 355 "zlamtsqr.f"
	    ztpmqrt_("R", "C", m, &kk, k, &c__0, nb, &a[ii + a_dim1], lda, &t[
		    (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, 
		    &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);
#line 358 "zlamtsqr.f"
	} else {
#line 359 "zlamtsqr.f"
	    ii = *n + 1;
#line 360 "zlamtsqr.f"
	}

#line 362 "zlamtsqr.f"
	i__1 = *mb + 1;
#line 362 "zlamtsqr.f"
	i__2 = -(*mb - *k);
#line 362 "zlamtsqr.f"
	for (i__ = ii - (*mb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__2) {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 366 "zlamtsqr.f"
	    --ctr;
#line 367 "zlamtsqr.f"
	    i__3 = *mb - *k;
#line 367 "zlamtsqr.f"
	    ztpmqrt_("R", "C", m, &i__3, k, &c__0, nb, &a[i__ + a_dim1], lda, 
		    &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 371 "zlamtsqr.f"
	}

/*         Multiply Q to the first block of C (1:M,1:MB) */

#line 375 "zlamtsqr.f"
	zgemqrt_("R", "C", m, mb, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 378 "zlamtsqr.f"
    } else if (right && notran) {

/*         Multiply Q to the first block of C */

#line 382 "zlamtsqr.f"
	kk = (*n - *k) % (*mb - *k);
#line 383 "zlamtsqr.f"
	ii = *n - kk + 1;
#line 384 "zlamtsqr.f"
	ctr = 1;
#line 385 "zlamtsqr.f"
	zgemqrt_("R", "N", m, mb, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], 
		ldt, &c__[c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		ftnlen)1);

#line 388 "zlamtsqr.f"
	i__2 = ii - *mb + *k;
#line 388 "zlamtsqr.f"
	i__1 = *mb - *k;
#line 388 "zlamtsqr.f"
	for (i__ = *mb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
		 {

/*         Multiply Q to the current block of C (1:M,I:I+MB) */

#line 392 "zlamtsqr.f"
	    i__3 = *mb - *k;
#line 392 "zlamtsqr.f"
	    ztpmqrt_("R", "N", m, &i__3, k, &c__0, nb, &a[i__ + a_dim1], lda, 
		    &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], 
		    ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info, (ftnlen)
		    1, (ftnlen)1);
#line 395 "zlamtsqr.f"
	    ++ctr;

#line 397 "zlamtsqr.f"
	}
#line 398 "zlamtsqr.f"
	if (ii <= *n) {

/*         Multiply Q to the last block of C */

#line 402 "zlamtsqr.f"
	    ztpmqrt_("R", "N", m, &kk, k, &c__0, nb, &a[ii + a_dim1], lda, &t[
		    (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, 
		    &c__[ii * c_dim1 + 1], ldc, &work[1], info, (ftnlen)1, (
		    ftnlen)1);

#line 406 "zlamtsqr.f"
	}

#line 408 "zlamtsqr.f"
    }

#line 410 "zlamtsqr.f"
    work[1].r = (doublereal) lw, work[1].i = 0.;
#line 411 "zlamtsqr.f"
    return 0;

/*     End of ZLAMTSQR */

} /* zlamtsqr_ */

