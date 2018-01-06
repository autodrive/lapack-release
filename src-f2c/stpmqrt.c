#line 1 "stpmqrt.f"
/* stpmqrt.f -- translated by f2c (version 20100827).
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

#line 1 "stpmqrt.f"
/* > \brief \b STPMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STPMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpmqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpmqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpmqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, */
/*                           A, LDA, B, LDB, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), */
/*      $          WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPMQRT applies a real orthogonal matrix Q obtained from a */
/* > "triangular-pentagonal" real block reflector H to a general */
/* > real matrix C, which consists of two blocks A and B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q^T from the Left; */
/* >          = 'R': apply Q or Q^T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q^T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The order of the trapezoidal part of V. */
/* >          K >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CTPQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is REAL array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CTPQRT in B.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDV >= max(1,M); */
/* >          if SIDE = 'R', LDV >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is REAL array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CTPQRT, stored as a NB-by-K matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension */
/* >          (LDA,N) if SIDE = 'L' or */
/* >          (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDC >= max(1,K); */
/* >          If SIDE = 'R', LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. */
/* >          LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array. The dimension of WORK is */
/* >           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'. */
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

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The columns of the pentagonal matrix V contain the elementary reflectors */
/* >  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a */
/* >  trapezoidal block V2: */
/* > */
/* >        V = [V1] */
/* >            [V2]. */
/* > */
/* >  The size of the trapezoidal block V2 is determined by the parameter L, */
/* >  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L */
/* >  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular; */
/* >  if L=0, there is no trapezoidal block, hence V = V1 is rectangular. */
/* > */
/* >  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K. */
/* >                      [B] */
/* > */
/* >  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K. */
/* > */
/* >  The real orthogonal matrix Q is formed from V and T. */
/* > */
/* >  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C. */
/* > */
/* >  If TRANS='T' and SIDE='L', C is on exit replaced with Q^T * C. */
/* > */
/* >  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q. */
/* > */
/* >  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q^T. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int stpmqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *l, integer *nb, doublereal *v, integer *ldv, 
	doublereal *t, integer *ldt, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *work, integer *info, ftnlen side_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, a_dim1, a_offset, b_dim1, b_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, lb, mb, kf, ldaq;
    static logical left, tran;
    static integer ldvq;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), stprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical notran;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     .. Test the input arguments .. */

#line 254 "stpmqrt.f"
    /* Parameter adjustments */
#line 254 "stpmqrt.f"
    v_dim1 = *ldv;
#line 254 "stpmqrt.f"
    v_offset = 1 + v_dim1;
#line 254 "stpmqrt.f"
    v -= v_offset;
#line 254 "stpmqrt.f"
    t_dim1 = *ldt;
#line 254 "stpmqrt.f"
    t_offset = 1 + t_dim1;
#line 254 "stpmqrt.f"
    t -= t_offset;
#line 254 "stpmqrt.f"
    a_dim1 = *lda;
#line 254 "stpmqrt.f"
    a_offset = 1 + a_dim1;
#line 254 "stpmqrt.f"
    a -= a_offset;
#line 254 "stpmqrt.f"
    b_dim1 = *ldb;
#line 254 "stpmqrt.f"
    b_offset = 1 + b_dim1;
#line 254 "stpmqrt.f"
    b -= b_offset;
#line 254 "stpmqrt.f"
    --work;
#line 254 "stpmqrt.f"

#line 254 "stpmqrt.f"
    /* Function Body */
#line 254 "stpmqrt.f"
    *info = 0;
#line 255 "stpmqrt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 256 "stpmqrt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 257 "stpmqrt.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 258 "stpmqrt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 260 "stpmqrt.f"
    if (left) {
#line 261 "stpmqrt.f"
	ldvq = max(1,*m);
#line 262 "stpmqrt.f"
	ldaq = max(1,*k);
#line 263 "stpmqrt.f"
    } else if (right) {
#line 264 "stpmqrt.f"
	ldvq = max(1,*n);
#line 265 "stpmqrt.f"
	ldaq = max(1,*m);
#line 266 "stpmqrt.f"
    }
#line 267 "stpmqrt.f"
    if (! left && ! right) {
#line 268 "stpmqrt.f"
	*info = -1;
#line 269 "stpmqrt.f"
    } else if (! tran && ! notran) {
#line 270 "stpmqrt.f"
	*info = -2;
#line 271 "stpmqrt.f"
    } else if (*m < 0) {
#line 272 "stpmqrt.f"
	*info = -3;
#line 273 "stpmqrt.f"
    } else if (*n < 0) {
#line 274 "stpmqrt.f"
	*info = -4;
#line 275 "stpmqrt.f"
    } else if (*k < 0) {
#line 276 "stpmqrt.f"
	*info = -5;
#line 277 "stpmqrt.f"
    } else if (*l < 0 || *l > *k) {
#line 278 "stpmqrt.f"
	*info = -6;
#line 279 "stpmqrt.f"
    } else if (*nb < 1 || *nb > *k && *k > 0) {
#line 280 "stpmqrt.f"
	*info = -7;
#line 281 "stpmqrt.f"
    } else if (*ldv < ldvq) {
#line 282 "stpmqrt.f"
	*info = -9;
#line 283 "stpmqrt.f"
    } else if (*ldt < *nb) {
#line 284 "stpmqrt.f"
	*info = -11;
#line 285 "stpmqrt.f"
    } else if (*lda < ldaq) {
#line 286 "stpmqrt.f"
	*info = -13;
#line 287 "stpmqrt.f"
    } else if (*ldb < max(1,*m)) {
#line 288 "stpmqrt.f"
	*info = -15;
#line 289 "stpmqrt.f"
    }

#line 291 "stpmqrt.f"
    if (*info != 0) {
#line 292 "stpmqrt.f"
	i__1 = -(*info);
#line 292 "stpmqrt.f"
	xerbla_("STPMQRT", &i__1, (ftnlen)7);
#line 293 "stpmqrt.f"
	return 0;
#line 294 "stpmqrt.f"
    }

/*     .. Quick return if possible .. */

#line 298 "stpmqrt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 298 "stpmqrt.f"
	return 0;
#line 298 "stpmqrt.f"
    }

#line 300 "stpmqrt.f"
    if (left && tran) {

#line 302 "stpmqrt.f"
	i__1 = *k;
#line 302 "stpmqrt.f"
	i__2 = *nb;
#line 302 "stpmqrt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 303 "stpmqrt.f"
	    i__3 = *nb, i__4 = *k - i__ + 1;
#line 303 "stpmqrt.f"
	    ib = min(i__3,i__4);
/* Computing MIN */
#line 304 "stpmqrt.f"
	    i__3 = *m - *l + i__ + ib - 1;
#line 304 "stpmqrt.f"
	    mb = min(i__3,*m);
#line 305 "stpmqrt.f"
	    if (i__ >= *l) {
#line 306 "stpmqrt.f"
		lb = 0;
#line 307 "stpmqrt.f"
	    } else {
#line 308 "stpmqrt.f"
		lb = mb - *m + *l - i__ + 1;
#line 309 "stpmqrt.f"
	    }
#line 310 "stpmqrt.f"
	    stprfb_("L", "T", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 313 "stpmqrt.f"
	}

#line 315 "stpmqrt.f"
    } else if (right && notran) {

#line 317 "stpmqrt.f"
	i__2 = *k;
#line 317 "stpmqrt.f"
	i__1 = *nb;
#line 317 "stpmqrt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 318 "stpmqrt.f"
	    i__3 = *nb, i__4 = *k - i__ + 1;
#line 318 "stpmqrt.f"
	    ib = min(i__3,i__4);
/* Computing MIN */
#line 319 "stpmqrt.f"
	    i__3 = *n - *l + i__ + ib - 1;
#line 319 "stpmqrt.f"
	    mb = min(i__3,*n);
#line 320 "stpmqrt.f"
	    if (i__ >= *l) {
#line 321 "stpmqrt.f"
		lb = 0;
#line 322 "stpmqrt.f"
	    } else {
#line 323 "stpmqrt.f"
		lb = mb - *n + *l - i__ + 1;
#line 324 "stpmqrt.f"
	    }
#line 325 "stpmqrt.f"
	    stprfb_("R", "N", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
#line 328 "stpmqrt.f"
	}

#line 330 "stpmqrt.f"
    } else if (left && notran) {

#line 332 "stpmqrt.f"
	kf = (*k - 1) / *nb * *nb + 1;
#line 333 "stpmqrt.f"
	i__1 = -(*nb);
#line 333 "stpmqrt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 334 "stpmqrt.f"
	    i__2 = *nb, i__3 = *k - i__ + 1;
#line 334 "stpmqrt.f"
	    ib = min(i__2,i__3);
/* Computing MIN */
#line 335 "stpmqrt.f"
	    i__2 = *m - *l + i__ + ib - 1;
#line 335 "stpmqrt.f"
	    mb = min(i__2,*m);
#line 336 "stpmqrt.f"
	    if (i__ >= *l) {
#line 337 "stpmqrt.f"
		lb = 0;
#line 338 "stpmqrt.f"
	    } else {
#line 339 "stpmqrt.f"
		lb = mb - *m + *l - i__ + 1;
#line 340 "stpmqrt.f"
	    }
#line 341 "stpmqrt.f"
	    stprfb_("L", "N", "F", "C", &mb, n, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &
		    b[b_offset], ldb, &work[1], &ib, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 344 "stpmqrt.f"
	}

#line 346 "stpmqrt.f"
    } else if (right && tran) {

#line 348 "stpmqrt.f"
	kf = (*k - 1) / *nb * *nb + 1;
#line 349 "stpmqrt.f"
	i__1 = -(*nb);
#line 349 "stpmqrt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 350 "stpmqrt.f"
	    i__2 = *nb, i__3 = *k - i__ + 1;
#line 350 "stpmqrt.f"
	    ib = min(i__2,i__3);
/* Computing MIN */
#line 351 "stpmqrt.f"
	    i__2 = *n - *l + i__ + ib - 1;
#line 351 "stpmqrt.f"
	    mb = min(i__2,*n);
#line 352 "stpmqrt.f"
	    if (i__ >= *l) {
#line 353 "stpmqrt.f"
		lb = 0;
#line 354 "stpmqrt.f"
	    } else {
#line 355 "stpmqrt.f"
		lb = mb - *n + *l - i__ + 1;
#line 356 "stpmqrt.f"
	    }
#line 357 "stpmqrt.f"
	    stprfb_("R", "T", "F", "C", m, &mb, &ib, &lb, &v[i__ * v_dim1 + 1]
		    , ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], 
		    lda, &b[b_offset], ldb, &work[1], m, (ftnlen)1, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);
#line 360 "stpmqrt.f"
	}

#line 362 "stpmqrt.f"
    }

#line 364 "stpmqrt.f"
    return 0;

/*     End of STPMQRT */

} /* stpmqrt_ */

