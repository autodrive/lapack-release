#line 1 "dgemlqt.f"
/* dgemlqt.f -- translated by f2c (version 20100827).
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

#line 1 "dgemlqt.f"
/* > \brief \b DGEMLQT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEMLQT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgemlqt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgemlqt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgemlqt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT, */
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
/* > DGEMLQT overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q C            C Q */
/* > TRANS = 'T':   Q**T C            C Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of K */
/* > elementary reflectors: */
/* > */
/* >       Q = H(1) H(2) . . . H(K) = I - V T V**T */
/* > */
/* > generated using the compact WY representation as returned by DGELQT. */
/* > */
/* > Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'. */
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
/* >          = 'C':  Transpose, apply Q**T. */
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
/* >          V is DOUBLE PRECISION array, dimension */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,K) */
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
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q. */
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
/* >          WORK is DOUBLE PRECISION array. The dimension of */
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
/* Subroutine */ int dgemlqt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *mb, doublereal *v, integer *ldv, doublereal *t, 
	integer *ldt, doublereal *c__, integer *ldc, doublereal *work, 
	integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, ib, kf;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
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

#line 205 "dgemlqt.f"
    /* Parameter adjustments */
#line 205 "dgemlqt.f"
    v_dim1 = *ldv;
#line 205 "dgemlqt.f"
    v_offset = 1 + v_dim1;
#line 205 "dgemlqt.f"
    v -= v_offset;
#line 205 "dgemlqt.f"
    t_dim1 = *ldt;
#line 205 "dgemlqt.f"
    t_offset = 1 + t_dim1;
#line 205 "dgemlqt.f"
    t -= t_offset;
#line 205 "dgemlqt.f"
    c_dim1 = *ldc;
#line 205 "dgemlqt.f"
    c_offset = 1 + c_dim1;
#line 205 "dgemlqt.f"
    c__ -= c_offset;
#line 205 "dgemlqt.f"
    --work;
#line 205 "dgemlqt.f"

#line 205 "dgemlqt.f"
    /* Function Body */
#line 205 "dgemlqt.f"
    *info = 0;
#line 206 "dgemlqt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 207 "dgemlqt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 208 "dgemlqt.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 209 "dgemlqt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 211 "dgemlqt.f"
    if (left) {
#line 212 "dgemlqt.f"
	ldwork = max(1,*n);
#line 213 "dgemlqt.f"
    } else if (right) {
#line 214 "dgemlqt.f"
	ldwork = max(1,*m);
#line 215 "dgemlqt.f"
    }
#line 216 "dgemlqt.f"
    if (! left && ! right) {
#line 217 "dgemlqt.f"
	*info = -1;
#line 218 "dgemlqt.f"
    } else if (! tran && ! notran) {
#line 219 "dgemlqt.f"
	*info = -2;
#line 220 "dgemlqt.f"
    } else if (*m < 0) {
#line 221 "dgemlqt.f"
	*info = -3;
#line 222 "dgemlqt.f"
    } else if (*n < 0) {
#line 223 "dgemlqt.f"
	*info = -4;
#line 224 "dgemlqt.f"
    } else if (*k < 0) {
#line 225 "dgemlqt.f"
	*info = -5;
#line 226 "dgemlqt.f"
    } else if (*mb < 1 || *mb > *k && *k > 0) {
#line 227 "dgemlqt.f"
	*info = -6;
#line 228 "dgemlqt.f"
    } else if (*ldv < max(1,*k)) {
#line 229 "dgemlqt.f"
	*info = -8;
#line 230 "dgemlqt.f"
    } else if (*ldt < *mb) {
#line 231 "dgemlqt.f"
	*info = -10;
#line 232 "dgemlqt.f"
    } else if (*ldc < max(1,*m)) {
#line 233 "dgemlqt.f"
	*info = -12;
#line 234 "dgemlqt.f"
    }

#line 236 "dgemlqt.f"
    if (*info != 0) {
#line 237 "dgemlqt.f"
	i__1 = -(*info);
#line 237 "dgemlqt.f"
	xerbla_("DGEMLQT", &i__1, (ftnlen)7);
#line 238 "dgemlqt.f"
	return 0;
#line 239 "dgemlqt.f"
    }

/*     .. Quick return if possible .. */

#line 243 "dgemlqt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 243 "dgemlqt.f"
	return 0;
#line 243 "dgemlqt.f"
    }

#line 245 "dgemlqt.f"
    if (left && notran) {

#line 247 "dgemlqt.f"
	i__1 = *k;
#line 247 "dgemlqt.f"
	i__2 = *mb;
#line 247 "dgemlqt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 248 "dgemlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 248 "dgemlqt.f"
	    ib = min(i__3,i__4);
#line 249 "dgemlqt.f"
	    i__3 = *m - i__ + 1;
#line 249 "dgemlqt.f"
	    dlarfb_("L", "T", "F", "R", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 252 "dgemlqt.f"
	}

#line 254 "dgemlqt.f"
    } else if (right && tran) {

#line 256 "dgemlqt.f"
	i__2 = *k;
#line 256 "dgemlqt.f"
	i__1 = *mb;
#line 256 "dgemlqt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 257 "dgemlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 257 "dgemlqt.f"
	    ib = min(i__3,i__4);
#line 258 "dgemlqt.f"
	    i__3 = *n - i__ + 1;
#line 258 "dgemlqt.f"
	    dlarfb_("R", "N", "F", "R", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 261 "dgemlqt.f"
	}

#line 263 "dgemlqt.f"
    } else if (left && tran) {

#line 265 "dgemlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 266 "dgemlqt.f"
	i__1 = -(*mb);
#line 266 "dgemlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 267 "dgemlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 267 "dgemlqt.f"
	    ib = min(i__2,i__3);
#line 268 "dgemlqt.f"
	    i__2 = *m - i__ + 1;
#line 268 "dgemlqt.f"
	    dlarfb_("L", "N", "F", "R", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 271 "dgemlqt.f"
	}

#line 273 "dgemlqt.f"
    } else if (right && notran) {

#line 275 "dgemlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 276 "dgemlqt.f"
	i__1 = -(*mb);
#line 276 "dgemlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 277 "dgemlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 277 "dgemlqt.f"
	    ib = min(i__2,i__3);
#line 278 "dgemlqt.f"
	    i__2 = *n - i__ + 1;
#line 278 "dgemlqt.f"
	    dlarfb_("R", "T", "F", "R", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 281 "dgemlqt.f"
	}

#line 283 "dgemlqt.f"
    }

#line 285 "dgemlqt.f"
    return 0;

/*     End of DGEMLQT */

} /* dgemlqt_ */

