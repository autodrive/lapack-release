#line 1 "ctpmlqt.f"
/* ctpmlqt.f -- translated by f2c (version 20100827).
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

#line 1 "ctpmlqt.f"
/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPMLQT( SIDE, TRANS, M, N, K, L, MB, V, LDV, T, LDT, */
/*                           A, LDA, B, LDB, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, MB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            V( LDV, * ), A( LDA, * ), B( LDB, * ), */
/*      $                   T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPMQRT applies a complex orthogonal matrix Q obtained from a */
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
/* >          V is COMPLEX array, dimension (LDA,K) */
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
/* >          T is COMPLEX array, dimension (LDT,K) */
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
/* >          A is COMPLEX array, dimension */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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
/* >          WORK is COMPLEX array. The dimension of WORK is */
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
/* Subroutine */ int ctpmlqt_(char *side, char *trans, integer *m, integer *n,
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), ctprfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical notran;


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

#line 237 "ctpmlqt.f"
    /* Parameter adjustments */
#line 237 "ctpmlqt.f"
    v_dim1 = *ldv;
#line 237 "ctpmlqt.f"
    v_offset = 1 + v_dim1;
#line 237 "ctpmlqt.f"
    v -= v_offset;
#line 237 "ctpmlqt.f"
    t_dim1 = *ldt;
#line 237 "ctpmlqt.f"
    t_offset = 1 + t_dim1;
#line 237 "ctpmlqt.f"
    t -= t_offset;
#line 237 "ctpmlqt.f"
    a_dim1 = *lda;
#line 237 "ctpmlqt.f"
    a_offset = 1 + a_dim1;
#line 237 "ctpmlqt.f"
    a -= a_offset;
#line 237 "ctpmlqt.f"
    b_dim1 = *ldb;
#line 237 "ctpmlqt.f"
    b_offset = 1 + b_dim1;
#line 237 "ctpmlqt.f"
    b -= b_offset;
#line 237 "ctpmlqt.f"
    --work;
#line 237 "ctpmlqt.f"

#line 237 "ctpmlqt.f"
    /* Function Body */
#line 237 "ctpmlqt.f"
    *info = 0;
#line 238 "ctpmlqt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 239 "ctpmlqt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 240 "ctpmlqt.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 241 "ctpmlqt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 243 "ctpmlqt.f"
    if (left) {
#line 244 "ctpmlqt.f"
	ldaq = max(1,*k);
#line 245 "ctpmlqt.f"
    } else if (right) {
#line 246 "ctpmlqt.f"
	ldaq = max(1,*m);
#line 247 "ctpmlqt.f"
    }
#line 248 "ctpmlqt.f"
    if (! left && ! right) {
#line 249 "ctpmlqt.f"
	*info = -1;
#line 250 "ctpmlqt.f"
    } else if (! tran && ! notran) {
#line 251 "ctpmlqt.f"
	*info = -2;
#line 252 "ctpmlqt.f"
    } else if (*m < 0) {
#line 253 "ctpmlqt.f"
	*info = -3;
#line 254 "ctpmlqt.f"
    } else if (*n < 0) {
#line 255 "ctpmlqt.f"
	*info = -4;
#line 256 "ctpmlqt.f"
    } else if (*k < 0) {
#line 257 "ctpmlqt.f"
	*info = -5;
#line 258 "ctpmlqt.f"
    } else if (*l < 0 || *l > *k) {
#line 259 "ctpmlqt.f"
	*info = -6;
#line 260 "ctpmlqt.f"
    } else if (*mb < 1 || *mb > *k && *k > 0) {
#line 261 "ctpmlqt.f"
	*info = -7;
#line 262 "ctpmlqt.f"
    } else if (*ldv < *k) {
#line 263 "ctpmlqt.f"
	*info = -9;
#line 264 "ctpmlqt.f"
    } else if (*ldt < *mb) {
#line 265 "ctpmlqt.f"
	*info = -11;
#line 266 "ctpmlqt.f"
    } else if (*lda < ldaq) {
#line 267 "ctpmlqt.f"
	*info = -13;
#line 268 "ctpmlqt.f"
    } else if (*ldb < max(1,*m)) {
#line 269 "ctpmlqt.f"
	*info = -15;
#line 270 "ctpmlqt.f"
    }

#line 272 "ctpmlqt.f"
    if (*info != 0) {
#line 273 "ctpmlqt.f"
	i__1 = -(*info);
#line 273 "ctpmlqt.f"
	xerbla_("CTPMLQT", &i__1, (ftnlen)7);
#line 274 "ctpmlqt.f"
	return 0;
#line 275 "ctpmlqt.f"
    }

/*     .. Quick return if possible .. */

#line 279 "ctpmlqt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 279 "ctpmlqt.f"
	return 0;
#line 279 "ctpmlqt.f"
    }

#line 281 "ctpmlqt.f"
    if (left && notran) {

#line 283 "ctpmlqt.f"
	i__1 = *k;
#line 283 "ctpmlqt.f"
	i__2 = *mb;
#line 283 "ctpmlqt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 284 "ctpmlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 284 "ctpmlqt.f"
	    ib = min(i__3,i__4);
/* Computing MIN */
#line 285 "ctpmlqt.f"
	    i__3 = *m - *l + i__ + ib - 1;
#line 285 "ctpmlqt.f"
	    nb = min(i__3,*m);
#line 286 "ctpmlqt.f"
	    if (i__ >= *l) {
#line 287 "ctpmlqt.f"
		lb = 0;
#line 288 "ctpmlqt.f"
	    } else {
#line 289 "ctpmlqt.f"
		lb = 0;
#line 290 "ctpmlqt.f"
	    }
#line 291 "ctpmlqt.f"
	    ctprfb_("L", "C", "F", "R", &nb, n, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &b[
		    b_offset], ldb, &work[1], &ib, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 294 "ctpmlqt.f"
	}

#line 296 "ctpmlqt.f"
    } else if (right && tran) {

#line 298 "ctpmlqt.f"
	i__2 = *k;
#line 298 "ctpmlqt.f"
	i__1 = *mb;
#line 298 "ctpmlqt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 299 "ctpmlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 299 "ctpmlqt.f"
	    ib = min(i__3,i__4);
/* Computing MIN */
#line 300 "ctpmlqt.f"
	    i__3 = *n - *l + i__ + ib - 1;
#line 300 "ctpmlqt.f"
	    nb = min(i__3,*n);
#line 301 "ctpmlqt.f"
	    if (i__ >= *l) {
#line 302 "ctpmlqt.f"
		lb = 0;
#line 303 "ctpmlqt.f"
	    } else {
#line 304 "ctpmlqt.f"
		lb = nb - *n + *l - i__ + 1;
#line 305 "ctpmlqt.f"
	    }
#line 306 "ctpmlqt.f"
	    ctprfb_("R", "N", "F", "R", m, &nb, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], lda,
		     &b[b_offset], ldb, &work[1], m, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 309 "ctpmlqt.f"
	}

#line 311 "ctpmlqt.f"
    } else if (left && tran) {

#line 313 "ctpmlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 314 "ctpmlqt.f"
	i__1 = -(*mb);
#line 314 "ctpmlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 315 "ctpmlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 315 "ctpmlqt.f"
	    ib = min(i__2,i__3);
/* Computing MIN */
#line 316 "ctpmlqt.f"
	    i__2 = *m - *l + i__ + ib - 1;
#line 316 "ctpmlqt.f"
	    nb = min(i__2,*m);
#line 317 "ctpmlqt.f"
	    if (i__ >= *l) {
#line 318 "ctpmlqt.f"
		lb = 0;
#line 319 "ctpmlqt.f"
	    } else {
#line 320 "ctpmlqt.f"
		lb = 0;
#line 321 "ctpmlqt.f"
	    }
#line 322 "ctpmlqt.f"
	    ctprfb_("L", "N", "F", "R", &nb, n, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &b[
		    b_offset], ldb, &work[1], &ib, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 325 "ctpmlqt.f"
	}

#line 327 "ctpmlqt.f"
    } else if (right && notran) {

#line 329 "ctpmlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 330 "ctpmlqt.f"
	i__1 = -(*mb);
#line 330 "ctpmlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 331 "ctpmlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 331 "ctpmlqt.f"
	    ib = min(i__2,i__3);
/* Computing MIN */
#line 332 "ctpmlqt.f"
	    i__2 = *n - *l + i__ + ib - 1;
#line 332 "ctpmlqt.f"
	    nb = min(i__2,*n);
#line 333 "ctpmlqt.f"
	    if (i__ >= *l) {
#line 334 "ctpmlqt.f"
		lb = 0;
#line 335 "ctpmlqt.f"
	    } else {
#line 336 "ctpmlqt.f"
		lb = nb - *n + *l - i__ + 1;
#line 337 "ctpmlqt.f"
	    }
#line 338 "ctpmlqt.f"
	    ctprfb_("R", "C", "F", "R", m, &nb, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], lda,
		     &b[b_offset], ldb, &work[1], m, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 341 "ctpmlqt.f"
	}

#line 343 "ctpmlqt.f"
    }

#line 345 "ctpmlqt.f"
    return 0;

/*     End of CTPMLQT */

} /* ctpmlqt_ */

