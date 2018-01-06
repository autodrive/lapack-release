#line 1 "ztpmlqt.f"
/* ztpmlqt.f -- translated by f2c (version 20100827).
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

#line 1 "ztpmlqt.f"
/* > \brief \b ZTPMLQT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpmlqt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpmlqt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpmlqt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPMLQT( SIDE, TRANS, M, N, K, L, MB, V, LDV, T, LDT, */
/*                           A, LDA, B, LDB, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, MB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         V( LDV, * ), A( LDA, * ), B( LDB, * ), */
/*      $                   T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPMQRT applies a complex orthogonal matrix Q obtained from a */
/* > "triangular-pentagonal" real block reflector H to a general */
/* > real matrix C, which consists of two blocks A and B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**C from the Left; */
/* >          = 'R': apply Q or Q**C from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Transpose, apply Q**C. */
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
/* > \param[in] MB */
/* > \verbatim */
/* >          MB is INTEGER */
/* >          The block size used for the storage of T.  K >= MB >= 1. */
/* >          This must be the same value of MB used to generate T */
/* >          in DTPLQT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (LDA,K) */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DTPLQT in B.  See Further Details. */
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
/* >          T is COMPLEX*16 array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by DTPLQT, stored as a MB-by-K matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= MB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension */
/* >          (LDA,N) if SIDE = 'L' or */
/* >          (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          Q*C or Q**C*C or C*Q or C*Q**C.  See Further Details. */
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
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          Q*C or Q**C*C or C*Q or C*Q**C.  See Further Details. */
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
/* >          WORK is COMPLEX*16 array. The dimension of WORK is */
/* >           N*MB if SIDE = 'L', or  M*MB if SIDE = 'R'. */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The columns of the pentagonal matrix V contain the elementary reflectors */
/* >  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a */
/* >  trapezoidal block V2: */
/* > */
/* >        V = [V1] [V2]. */
/* > */
/* > */
/* >  The size of the trapezoidal block V2 is determined by the parameter L, */
/* >  where 0 <= L <= K; V2 is lower trapezoidal, consisting of the first L */
/* >  rows of a K-by-K upper triangular matrix.  If L=K, V2 is lower triangular; */
/* >  if L=0, there is no trapezoidal block, hence V = V1 is rectangular. */
/* > */
/* >  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is K-by-M. */
/* >                      [B] */
/* > */
/* >  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is K-by-N. */
/* > */
/* >  The real orthogonal matrix Q is formed from V and T. */
/* > */
/* >  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C. */
/* > */
/* >  If TRANS='C' and SIDE='L', C is on exit replaced with Q**C * C. */
/* > */
/* >  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q. */
/* > */
/* >  If TRANS='C' and SIDE='R', C is on exit replaced with C * Q**C. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztpmlqt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *l, integer *mb, doublecomplex *v, integer *ldv, 
	doublecomplex *t, integer *ldt, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *work, integer *info, 
	ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, a_dim1, a_offset, b_dim1, b_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, lb, nb, kf, ldaq;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    extern /* Subroutine */ int ztprfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 254 "ztpmlqt.f"
    /* Parameter adjustments */
#line 254 "ztpmlqt.f"
    v_dim1 = *ldv;
#line 254 "ztpmlqt.f"
    v_offset = 1 + v_dim1;
#line 254 "ztpmlqt.f"
    v -= v_offset;
#line 254 "ztpmlqt.f"
    t_dim1 = *ldt;
#line 254 "ztpmlqt.f"
    t_offset = 1 + t_dim1;
#line 254 "ztpmlqt.f"
    t -= t_offset;
#line 254 "ztpmlqt.f"
    a_dim1 = *lda;
#line 254 "ztpmlqt.f"
    a_offset = 1 + a_dim1;
#line 254 "ztpmlqt.f"
    a -= a_offset;
#line 254 "ztpmlqt.f"
    b_dim1 = *ldb;
#line 254 "ztpmlqt.f"
    b_offset = 1 + b_dim1;
#line 254 "ztpmlqt.f"
    b -= b_offset;
#line 254 "ztpmlqt.f"
    --work;
#line 254 "ztpmlqt.f"

#line 254 "ztpmlqt.f"
    /* Function Body */
#line 254 "ztpmlqt.f"
    *info = 0;
#line 255 "ztpmlqt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 256 "ztpmlqt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 257 "ztpmlqt.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 258 "ztpmlqt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 260 "ztpmlqt.f"
    if (left) {
#line 261 "ztpmlqt.f"
	ldaq = max(1,*k);
#line 262 "ztpmlqt.f"
    } else if (right) {
#line 263 "ztpmlqt.f"
	ldaq = max(1,*m);
#line 264 "ztpmlqt.f"
    }
#line 265 "ztpmlqt.f"
    if (! left && ! right) {
#line 266 "ztpmlqt.f"
	*info = -1;
#line 267 "ztpmlqt.f"
    } else if (! tran && ! notran) {
#line 268 "ztpmlqt.f"
	*info = -2;
#line 269 "ztpmlqt.f"
    } else if (*m < 0) {
#line 270 "ztpmlqt.f"
	*info = -3;
#line 271 "ztpmlqt.f"
    } else if (*n < 0) {
#line 272 "ztpmlqt.f"
	*info = -4;
#line 273 "ztpmlqt.f"
    } else if (*k < 0) {
#line 274 "ztpmlqt.f"
	*info = -5;
#line 275 "ztpmlqt.f"
    } else if (*l < 0 || *l > *k) {
#line 276 "ztpmlqt.f"
	*info = -6;
#line 277 "ztpmlqt.f"
    } else if (*mb < 1 || *mb > *k && *k > 0) {
#line 278 "ztpmlqt.f"
	*info = -7;
#line 279 "ztpmlqt.f"
    } else if (*ldv < *k) {
#line 280 "ztpmlqt.f"
	*info = -9;
#line 281 "ztpmlqt.f"
    } else if (*ldt < *mb) {
#line 282 "ztpmlqt.f"
	*info = -11;
#line 283 "ztpmlqt.f"
    } else if (*lda < ldaq) {
#line 284 "ztpmlqt.f"
	*info = -13;
#line 285 "ztpmlqt.f"
    } else if (*ldb < max(1,*m)) {
#line 286 "ztpmlqt.f"
	*info = -15;
#line 287 "ztpmlqt.f"
    }

#line 289 "ztpmlqt.f"
    if (*info != 0) {
#line 290 "ztpmlqt.f"
	i__1 = -(*info);
#line 290 "ztpmlqt.f"
	xerbla_("ZTPMLQT", &i__1, (ftnlen)7);
#line 291 "ztpmlqt.f"
	return 0;
#line 292 "ztpmlqt.f"
    }

/*     .. Quick return if possible .. */

#line 296 "ztpmlqt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 296 "ztpmlqt.f"
	return 0;
#line 296 "ztpmlqt.f"
    }

#line 298 "ztpmlqt.f"
    if (left && notran) {

#line 300 "ztpmlqt.f"
	i__1 = *k;
#line 300 "ztpmlqt.f"
	i__2 = *mb;
#line 300 "ztpmlqt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 301 "ztpmlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 301 "ztpmlqt.f"
	    ib = min(i__3,i__4);
/* Computing MIN */
#line 302 "ztpmlqt.f"
	    i__3 = *m - *l + i__ + ib - 1;
#line 302 "ztpmlqt.f"
	    nb = min(i__3,*m);
#line 303 "ztpmlqt.f"
	    if (i__ >= *l) {
#line 304 "ztpmlqt.f"
		lb = 0;
#line 305 "ztpmlqt.f"
	    } else {
#line 306 "ztpmlqt.f"
		lb = 0;
#line 307 "ztpmlqt.f"
	    }
#line 308 "ztpmlqt.f"
	    ztprfb_("L", "C", "F", "R", &nb, n, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &b[
		    b_offset], ldb, &work[1], &ib, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 311 "ztpmlqt.f"
	}

#line 313 "ztpmlqt.f"
    } else if (right && tran) {

#line 315 "ztpmlqt.f"
	i__2 = *k;
#line 315 "ztpmlqt.f"
	i__1 = *mb;
#line 315 "ztpmlqt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 316 "ztpmlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 316 "ztpmlqt.f"
	    ib = min(i__3,i__4);
/* Computing MIN */
#line 317 "ztpmlqt.f"
	    i__3 = *n - *l + i__ + ib - 1;
#line 317 "ztpmlqt.f"
	    nb = min(i__3,*n);
#line 318 "ztpmlqt.f"
	    if (i__ >= *l) {
#line 319 "ztpmlqt.f"
		lb = 0;
#line 320 "ztpmlqt.f"
	    } else {
#line 321 "ztpmlqt.f"
		lb = nb - *n + *l - i__ + 1;
#line 322 "ztpmlqt.f"
	    }
#line 323 "ztpmlqt.f"
	    ztprfb_("R", "N", "F", "R", m, &nb, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], lda,
		     &b[b_offset], ldb, &work[1], m, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 326 "ztpmlqt.f"
	}

#line 328 "ztpmlqt.f"
    } else if (left && tran) {

#line 330 "ztpmlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 331 "ztpmlqt.f"
	i__1 = -(*mb);
#line 331 "ztpmlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 332 "ztpmlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 332 "ztpmlqt.f"
	    ib = min(i__2,i__3);
/* Computing MIN */
#line 333 "ztpmlqt.f"
	    i__2 = *m - *l + i__ + ib - 1;
#line 333 "ztpmlqt.f"
	    nb = min(i__2,*m);
#line 334 "ztpmlqt.f"
	    if (i__ >= *l) {
#line 335 "ztpmlqt.f"
		lb = 0;
#line 336 "ztpmlqt.f"
	    } else {
#line 337 "ztpmlqt.f"
		lb = 0;
#line 338 "ztpmlqt.f"
	    }
#line 339 "ztpmlqt.f"
	    ztprfb_("L", "N", "F", "R", &nb, n, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &b[
		    b_offset], ldb, &work[1], &ib, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 342 "ztpmlqt.f"
	}

#line 344 "ztpmlqt.f"
    } else if (right && notran) {

#line 346 "ztpmlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 347 "ztpmlqt.f"
	i__1 = -(*mb);
#line 347 "ztpmlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 348 "ztpmlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 348 "ztpmlqt.f"
	    ib = min(i__2,i__3);
/* Computing MIN */
#line 349 "ztpmlqt.f"
	    i__2 = *n - *l + i__ + ib - 1;
#line 349 "ztpmlqt.f"
	    nb = min(i__2,*n);
#line 350 "ztpmlqt.f"
	    if (i__ >= *l) {
#line 351 "ztpmlqt.f"
		lb = 0;
#line 352 "ztpmlqt.f"
	    } else {
#line 353 "ztpmlqt.f"
		lb = nb - *n + *l - i__ + 1;
#line 354 "ztpmlqt.f"
	    }
#line 355 "ztpmlqt.f"
	    ztprfb_("R", "C", "F", "R", m, &nb, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], lda,
		     &b[b_offset], ldb, &work[1], m, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 358 "ztpmlqt.f"
	}

#line 360 "ztpmlqt.f"
    }

#line 362 "ztpmlqt.f"
    return 0;

/*     End of ZTPMLQT */

} /* ztpmlqt_ */

