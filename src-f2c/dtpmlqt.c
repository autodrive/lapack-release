#line 1 "dtpmlqt.f"
/* dtpmlqt.f -- translated by f2c (version 20100827).
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

#line 1 "dtpmlqt.f"
/* > \brief \b DTPMLQT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpmlqt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpmlqt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpmlqt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPMLQT( SIDE, TRANS, M, N, K, L, MB, V, LDV, T, LDT, */
/*                           A, LDA, B, LDB, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, MB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   V( LDV, * ), A( LDA, * ), B( LDB, * ), */
/*      $                   T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPMQRT applies a real orthogonal matrix Q obtained from a */
/* > "triangular-pentagonal" real block reflector H to a general */
/* > real matrix C, which consists of two blocks A and B. */
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
/* >          V is DOUBLE PRECISION array, dimension (LDA,K) */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,K) */
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
/* >          A is DOUBLE PRECISION array, dimension */
/* >          (LDA,N) if SIDE = 'L' or */
/* >          (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          Q*C or Q**T*C or C*Q or C*Q**T.  See Further Details. */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          Q*C or Q**T*C or C*Q or C*Q**T.  See Further Details. */
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
/* >          WORK is DOUBLE PRECISION array. The dimension of WORK is */
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

/* > \date November 2017 */

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
/* >  If TRANS='T' and SIDE='L', C is on exit replaced with Q**T * C. */
/* > */
/* >  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q. */
/* > */
/* >  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q**T. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtpmlqt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *l, integer *mb, doublereal *v, integer *ldv, 
	doublereal *t, integer *ldt, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *work, integer *info, ftnlen side_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, a_dim1, a_offset, b_dim1, b_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, lb, nb, kf, ldaq;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dtprfb_(
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

#line 254 "dtpmlqt.f"
    /* Parameter adjustments */
#line 254 "dtpmlqt.f"
    v_dim1 = *ldv;
#line 254 "dtpmlqt.f"
    v_offset = 1 + v_dim1;
#line 254 "dtpmlqt.f"
    v -= v_offset;
#line 254 "dtpmlqt.f"
    t_dim1 = *ldt;
#line 254 "dtpmlqt.f"
    t_offset = 1 + t_dim1;
#line 254 "dtpmlqt.f"
    t -= t_offset;
#line 254 "dtpmlqt.f"
    a_dim1 = *lda;
#line 254 "dtpmlqt.f"
    a_offset = 1 + a_dim1;
#line 254 "dtpmlqt.f"
    a -= a_offset;
#line 254 "dtpmlqt.f"
    b_dim1 = *ldb;
#line 254 "dtpmlqt.f"
    b_offset = 1 + b_dim1;
#line 254 "dtpmlqt.f"
    b -= b_offset;
#line 254 "dtpmlqt.f"
    --work;
#line 254 "dtpmlqt.f"

#line 254 "dtpmlqt.f"
    /* Function Body */
#line 254 "dtpmlqt.f"
    *info = 0;
#line 255 "dtpmlqt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 256 "dtpmlqt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 257 "dtpmlqt.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 258 "dtpmlqt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 260 "dtpmlqt.f"
    if (left) {
#line 261 "dtpmlqt.f"
	ldaq = max(1,*k);
#line 262 "dtpmlqt.f"
    } else if (right) {
#line 263 "dtpmlqt.f"
	ldaq = max(1,*m);
#line 264 "dtpmlqt.f"
    }
#line 265 "dtpmlqt.f"
    if (! left && ! right) {
#line 266 "dtpmlqt.f"
	*info = -1;
#line 267 "dtpmlqt.f"
    } else if (! tran && ! notran) {
#line 268 "dtpmlqt.f"
	*info = -2;
#line 269 "dtpmlqt.f"
    } else if (*m < 0) {
#line 270 "dtpmlqt.f"
	*info = -3;
#line 271 "dtpmlqt.f"
    } else if (*n < 0) {
#line 272 "dtpmlqt.f"
	*info = -4;
#line 273 "dtpmlqt.f"
    } else if (*k < 0) {
#line 274 "dtpmlqt.f"
	*info = -5;
#line 275 "dtpmlqt.f"
    } else if (*l < 0 || *l > *k) {
#line 276 "dtpmlqt.f"
	*info = -6;
#line 277 "dtpmlqt.f"
    } else if (*mb < 1 || *mb > *k && *k > 0) {
#line 278 "dtpmlqt.f"
	*info = -7;
#line 279 "dtpmlqt.f"
    } else if (*ldv < *k) {
#line 280 "dtpmlqt.f"
	*info = -9;
#line 281 "dtpmlqt.f"
    } else if (*ldt < *mb) {
#line 282 "dtpmlqt.f"
	*info = -11;
#line 283 "dtpmlqt.f"
    } else if (*lda < ldaq) {
#line 284 "dtpmlqt.f"
	*info = -13;
#line 285 "dtpmlqt.f"
    } else if (*ldb < max(1,*m)) {
#line 286 "dtpmlqt.f"
	*info = -15;
#line 287 "dtpmlqt.f"
    }

#line 289 "dtpmlqt.f"
    if (*info != 0) {
#line 290 "dtpmlqt.f"
	i__1 = -(*info);
#line 290 "dtpmlqt.f"
	xerbla_("DTPMLQT", &i__1, (ftnlen)7);
#line 291 "dtpmlqt.f"
	return 0;
#line 292 "dtpmlqt.f"
    }

/*     .. Quick return if possible .. */

#line 296 "dtpmlqt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 296 "dtpmlqt.f"
	return 0;
#line 296 "dtpmlqt.f"
    }

#line 298 "dtpmlqt.f"
    if (left && notran) {

#line 300 "dtpmlqt.f"
	i__1 = *k;
#line 300 "dtpmlqt.f"
	i__2 = *mb;
#line 300 "dtpmlqt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 301 "dtpmlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 301 "dtpmlqt.f"
	    ib = min(i__3,i__4);
/* Computing MIN */
#line 302 "dtpmlqt.f"
	    i__3 = *m - *l + i__ + ib - 1;
#line 302 "dtpmlqt.f"
	    nb = min(i__3,*m);
#line 303 "dtpmlqt.f"
	    if (i__ >= *l) {
#line 304 "dtpmlqt.f"
		lb = 0;
#line 305 "dtpmlqt.f"
	    } else {
#line 306 "dtpmlqt.f"
		lb = 0;
#line 307 "dtpmlqt.f"
	    }
#line 308 "dtpmlqt.f"
	    dtprfb_("L", "T", "F", "R", &nb, n, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &b[
		    b_offset], ldb, &work[1], &ib, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 311 "dtpmlqt.f"
	}

#line 313 "dtpmlqt.f"
    } else if (right && tran) {

#line 315 "dtpmlqt.f"
	i__2 = *k;
#line 315 "dtpmlqt.f"
	i__1 = *mb;
#line 315 "dtpmlqt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 316 "dtpmlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 316 "dtpmlqt.f"
	    ib = min(i__3,i__4);
/* Computing MIN */
#line 317 "dtpmlqt.f"
	    i__3 = *n - *l + i__ + ib - 1;
#line 317 "dtpmlqt.f"
	    nb = min(i__3,*n);
#line 318 "dtpmlqt.f"
	    if (i__ >= *l) {
#line 319 "dtpmlqt.f"
		lb = 0;
#line 320 "dtpmlqt.f"
	    } else {
#line 321 "dtpmlqt.f"
		lb = nb - *n + *l - i__ + 1;
#line 322 "dtpmlqt.f"
	    }
#line 323 "dtpmlqt.f"
	    dtprfb_("R", "N", "F", "R", m, &nb, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], lda,
		     &b[b_offset], ldb, &work[1], m, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 326 "dtpmlqt.f"
	}

#line 328 "dtpmlqt.f"
    } else if (left && tran) {

#line 330 "dtpmlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 331 "dtpmlqt.f"
	i__1 = -(*mb);
#line 331 "dtpmlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 332 "dtpmlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 332 "dtpmlqt.f"
	    ib = min(i__2,i__3);
/* Computing MIN */
#line 333 "dtpmlqt.f"
	    i__2 = *m - *l + i__ + ib - 1;
#line 333 "dtpmlqt.f"
	    nb = min(i__2,*m);
#line 334 "dtpmlqt.f"
	    if (i__ >= *l) {
#line 335 "dtpmlqt.f"
		lb = 0;
#line 336 "dtpmlqt.f"
	    } else {
#line 337 "dtpmlqt.f"
		lb = 0;
#line 338 "dtpmlqt.f"
	    }
#line 339 "dtpmlqt.f"
	    dtprfb_("L", "N", "F", "R", &nb, n, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &b[
		    b_offset], ldb, &work[1], &ib, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 342 "dtpmlqt.f"
	}

#line 344 "dtpmlqt.f"
    } else if (right && notran) {

#line 346 "dtpmlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 347 "dtpmlqt.f"
	i__1 = -(*mb);
#line 347 "dtpmlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 348 "dtpmlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 348 "dtpmlqt.f"
	    ib = min(i__2,i__3);
/* Computing MIN */
#line 349 "dtpmlqt.f"
	    i__2 = *n - *l + i__ + ib - 1;
#line 349 "dtpmlqt.f"
	    nb = min(i__2,*n);
#line 350 "dtpmlqt.f"
	    if (i__ >= *l) {
#line 351 "dtpmlqt.f"
		lb = 0;
#line 352 "dtpmlqt.f"
	    } else {
#line 353 "dtpmlqt.f"
		lb = nb - *n + *l - i__ + 1;
#line 354 "dtpmlqt.f"
	    }
#line 355 "dtpmlqt.f"
	    dtprfb_("R", "T", "F", "R", m, &nb, &ib, &lb, &v[i__ + v_dim1], 
		    ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], lda,
		     &b[b_offset], ldb, &work[1], m, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 358 "dtpmlqt.f"
	}

#line 360 "dtpmlqt.f"
    }

#line 362 "dtpmlqt.f"
    return 0;

/*     End of DTPMLQT */

} /* dtpmlqt_ */

