#line 1 "zgemlqt.f"
/* zgemlqt.f -- translated by f2c (version 20100827).
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

#line 1 "zgemlqt.f"
/* > \brief \b ZGEMLQT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEMLQT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgemlqt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgemlqt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgemlqt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT, */
/*                          C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDC, M, N, MB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEMLQT overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q C            C Q */
/* > TRANS = 'C':   Q**H C            C Q**H */
/* > */
/* > where Q is a complex orthogonal matrix defined as the product of K */
/* > elementary reflectors: */
/* > */
/* >       Q = H(1) H(2) . . . H(K) = I - V T V**H */
/* > */
/* > generated using the compact WY representation as returned by ZGELQT. */
/* > */
/* > Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'. */
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
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          If SIDE = 'L', M >= K >= 0; */
/* >          if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* >          MB is INTEGER */
/* >          The block size used for the storage of T.  K >= MB >= 1. */
/* >          This must be the same value of MB used to generate T */
/* >          in DGELQT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension */
/* >                               (LDV,M) if SIDE = 'L', */
/* >                               (LDV,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DGELQT in the first K rows of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. LDV >= max(1,K). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by DGELQT, stored as a MB-by-K matrix. */
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
/* >          On exit, C is overwritten by Q C, Q**H C, C Q**H or C Q. */
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
/* >          WORK is COMPLEX*16 array. The dimension of */
/* >          WORK is N*MB if SIDE = 'L', or  M*MB if SIDE = 'R'. */
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

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgemlqt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *mb, doublecomplex *v, integer *ldv, 
	doublecomplex *t, integer *ldt, doublecomplex *c__, integer *ldc, 
	doublecomplex *work, integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, ib, kf;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlarfb_(
	    char *, char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static logical notran;
    static integer ldwork;


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

#line 205 "zgemlqt.f"
    /* Parameter adjustments */
#line 205 "zgemlqt.f"
    v_dim1 = *ldv;
#line 205 "zgemlqt.f"
    v_offset = 1 + v_dim1;
#line 205 "zgemlqt.f"
    v -= v_offset;
#line 205 "zgemlqt.f"
    t_dim1 = *ldt;
#line 205 "zgemlqt.f"
    t_offset = 1 + t_dim1;
#line 205 "zgemlqt.f"
    t -= t_offset;
#line 205 "zgemlqt.f"
    c_dim1 = *ldc;
#line 205 "zgemlqt.f"
    c_offset = 1 + c_dim1;
#line 205 "zgemlqt.f"
    c__ -= c_offset;
#line 205 "zgemlqt.f"
    --work;
#line 205 "zgemlqt.f"

#line 205 "zgemlqt.f"
    /* Function Body */
#line 205 "zgemlqt.f"
    *info = 0;
#line 206 "zgemlqt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 207 "zgemlqt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 208 "zgemlqt.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 209 "zgemlqt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 211 "zgemlqt.f"
    if (left) {
#line 212 "zgemlqt.f"
	ldwork = max(1,*n);
#line 213 "zgemlqt.f"
    } else if (right) {
#line 214 "zgemlqt.f"
	ldwork = max(1,*m);
#line 215 "zgemlqt.f"
    }
#line 216 "zgemlqt.f"
    if (! left && ! right) {
#line 217 "zgemlqt.f"
	*info = -1;
#line 218 "zgemlqt.f"
    } else if (! tran && ! notran) {
#line 219 "zgemlqt.f"
	*info = -2;
#line 220 "zgemlqt.f"
    } else if (*m < 0) {
#line 221 "zgemlqt.f"
	*info = -3;
#line 222 "zgemlqt.f"
    } else if (*n < 0) {
#line 223 "zgemlqt.f"
	*info = -4;
#line 224 "zgemlqt.f"
    } else if (*k < 0) {
#line 225 "zgemlqt.f"
	*info = -5;
#line 226 "zgemlqt.f"
    } else if (*mb < 1 || *mb > *k && *k > 0) {
#line 227 "zgemlqt.f"
	*info = -6;
#line 228 "zgemlqt.f"
    } else if (*ldv < max(1,*k)) {
#line 229 "zgemlqt.f"
	*info = -8;
#line 230 "zgemlqt.f"
    } else if (*ldt < *mb) {
#line 231 "zgemlqt.f"
	*info = -10;
#line 232 "zgemlqt.f"
    } else if (*ldc < max(1,*m)) {
#line 233 "zgemlqt.f"
	*info = -12;
#line 234 "zgemlqt.f"
    }

#line 236 "zgemlqt.f"
    if (*info != 0) {
#line 237 "zgemlqt.f"
	i__1 = -(*info);
#line 237 "zgemlqt.f"
	xerbla_("ZGEMLQT", &i__1, (ftnlen)7);
#line 238 "zgemlqt.f"
	return 0;
#line 239 "zgemlqt.f"
    }

/*     .. Quick return if possible .. */

#line 243 "zgemlqt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 243 "zgemlqt.f"
	return 0;
#line 243 "zgemlqt.f"
    }

#line 245 "zgemlqt.f"
    if (left && notran) {

#line 247 "zgemlqt.f"
	i__1 = *k;
#line 247 "zgemlqt.f"
	i__2 = *mb;
#line 247 "zgemlqt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 248 "zgemlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 248 "zgemlqt.f"
	    ib = min(i__3,i__4);
#line 249 "zgemlqt.f"
	    i__3 = *m - i__ + 1;
#line 249 "zgemlqt.f"
	    zlarfb_("L", "C", "F", "R", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 252 "zgemlqt.f"
	}

#line 254 "zgemlqt.f"
    } else if (right && tran) {

#line 256 "zgemlqt.f"
	i__2 = *k;
#line 256 "zgemlqt.f"
	i__1 = *mb;
#line 256 "zgemlqt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 257 "zgemlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 257 "zgemlqt.f"
	    ib = min(i__3,i__4);
#line 258 "zgemlqt.f"
	    i__3 = *n - i__ + 1;
#line 258 "zgemlqt.f"
	    zlarfb_("R", "N", "F", "R", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 261 "zgemlqt.f"
	}

#line 263 "zgemlqt.f"
    } else if (left && tran) {

#line 265 "zgemlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 266 "zgemlqt.f"
	i__1 = -(*mb);
#line 266 "zgemlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 267 "zgemlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 267 "zgemlqt.f"
	    ib = min(i__2,i__3);
#line 268 "zgemlqt.f"
	    i__2 = *m - i__ + 1;
#line 268 "zgemlqt.f"
	    zlarfb_("L", "N", "F", "R", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 271 "zgemlqt.f"
	}

#line 273 "zgemlqt.f"
    } else if (right && notran) {

#line 275 "zgemlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 276 "zgemlqt.f"
	i__1 = -(*mb);
#line 276 "zgemlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 277 "zgemlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 277 "zgemlqt.f"
	    ib = min(i__2,i__3);
#line 278 "zgemlqt.f"
	    i__2 = *n - i__ + 1;
#line 278 "zgemlqt.f"
	    zlarfb_("R", "C", "F", "R", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 281 "zgemlqt.f"
	}

#line 283 "zgemlqt.f"
    }

#line 285 "zgemlqt.f"
    return 0;

/*     End of ZGEMLQT */

} /* zgemlqt_ */

