#line 1 "sgemlqt.f"
/* sgemlqt.f -- translated by f2c (version 20100827).
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

#line 1 "sgemlqt.f"
/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT, */
/*                          C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDC, M, N, MB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL      V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
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
/* >          V is REAL array, dimension */
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
/* >          T is REAL array, dimension (LDT,K) */
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
/* >          C is REAL array, dimension (LDC,N) */
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
/* >          WORK is REAL array. The dimension of */
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
/* Subroutine */ int sgemlqt_(char *side, char *trans, integer *m, integer *n,
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
    extern /* Subroutine */ int slarfb_(char *, char *, char *, char *, 
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

#line 188 "sgemlqt.f"
    /* Parameter adjustments */
#line 188 "sgemlqt.f"
    v_dim1 = *ldv;
#line 188 "sgemlqt.f"
    v_offset = 1 + v_dim1;
#line 188 "sgemlqt.f"
    v -= v_offset;
#line 188 "sgemlqt.f"
    t_dim1 = *ldt;
#line 188 "sgemlqt.f"
    t_offset = 1 + t_dim1;
#line 188 "sgemlqt.f"
    t -= t_offset;
#line 188 "sgemlqt.f"
    c_dim1 = *ldc;
#line 188 "sgemlqt.f"
    c_offset = 1 + c_dim1;
#line 188 "sgemlqt.f"
    c__ -= c_offset;
#line 188 "sgemlqt.f"
    --work;
#line 188 "sgemlqt.f"

#line 188 "sgemlqt.f"
    /* Function Body */
#line 188 "sgemlqt.f"
    *info = 0;
#line 189 "sgemlqt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 190 "sgemlqt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 191 "sgemlqt.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 192 "sgemlqt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 194 "sgemlqt.f"
    if (left) {
#line 195 "sgemlqt.f"
	ldwork = max(1,*n);
#line 196 "sgemlqt.f"
    } else if (right) {
#line 197 "sgemlqt.f"
	ldwork = max(1,*m);
#line 198 "sgemlqt.f"
    }
#line 199 "sgemlqt.f"
    if (! left && ! right) {
#line 200 "sgemlqt.f"
	*info = -1;
#line 201 "sgemlqt.f"
    } else if (! tran && ! notran) {
#line 202 "sgemlqt.f"
	*info = -2;
#line 203 "sgemlqt.f"
    } else if (*m < 0) {
#line 204 "sgemlqt.f"
	*info = -3;
#line 205 "sgemlqt.f"
    } else if (*n < 0) {
#line 206 "sgemlqt.f"
	*info = -4;
#line 207 "sgemlqt.f"
    } else if (*k < 0) {
#line 208 "sgemlqt.f"
	*info = -5;
#line 209 "sgemlqt.f"
    } else if (*mb < 1 || *mb > *k && *k > 0) {
#line 210 "sgemlqt.f"
	*info = -6;
#line 211 "sgemlqt.f"
    } else if (*ldv < max(1,*k)) {
#line 212 "sgemlqt.f"
	*info = -8;
#line 213 "sgemlqt.f"
    } else if (*ldt < *mb) {
#line 214 "sgemlqt.f"
	*info = -10;
#line 215 "sgemlqt.f"
    } else if (*ldc < max(1,*m)) {
#line 216 "sgemlqt.f"
	*info = -12;
#line 217 "sgemlqt.f"
    }

#line 219 "sgemlqt.f"
    if (*info != 0) {
#line 220 "sgemlqt.f"
	i__1 = -(*info);
#line 220 "sgemlqt.f"
	xerbla_("SGEMLQT", &i__1, (ftnlen)7);
#line 221 "sgemlqt.f"
	return 0;
#line 222 "sgemlqt.f"
    }

/*     .. Quick return if possible .. */

#line 226 "sgemlqt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 226 "sgemlqt.f"
	return 0;
#line 226 "sgemlqt.f"
    }

#line 228 "sgemlqt.f"
    if (left && notran) {

#line 230 "sgemlqt.f"
	i__1 = *k;
#line 230 "sgemlqt.f"
	i__2 = *mb;
#line 230 "sgemlqt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 231 "sgemlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 231 "sgemlqt.f"
	    ib = min(i__3,i__4);
#line 232 "sgemlqt.f"
	    i__3 = *m - i__ + 1;
#line 232 "sgemlqt.f"
	    slarfb_("L", "T", "F", "R", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 235 "sgemlqt.f"
	}

#line 237 "sgemlqt.f"
    } else if (right && tran) {

#line 239 "sgemlqt.f"
	i__2 = *k;
#line 239 "sgemlqt.f"
	i__1 = *mb;
#line 239 "sgemlqt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 240 "sgemlqt.f"
	    i__3 = *mb, i__4 = *k - i__ + 1;
#line 240 "sgemlqt.f"
	    ib = min(i__3,i__4);
#line 241 "sgemlqt.f"
	    i__3 = *n - i__ + 1;
#line 241 "sgemlqt.f"
	    slarfb_("R", "N", "F", "R", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 244 "sgemlqt.f"
	}

#line 246 "sgemlqt.f"
    } else if (left && tran) {

#line 248 "sgemlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 249 "sgemlqt.f"
	i__1 = -(*mb);
#line 249 "sgemlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 250 "sgemlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 250 "sgemlqt.f"
	    ib = min(i__2,i__3);
#line 251 "sgemlqt.f"
	    i__2 = *m - i__ + 1;
#line 251 "sgemlqt.f"
	    slarfb_("L", "N", "F", "R", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 254 "sgemlqt.f"
	}

#line 256 "sgemlqt.f"
    } else if (right && notran) {

#line 258 "sgemlqt.f"
	kf = (*k - 1) / *mb * *mb + 1;
#line 259 "sgemlqt.f"
	i__1 = -(*mb);
#line 259 "sgemlqt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 260 "sgemlqt.f"
	    i__2 = *mb, i__3 = *k - i__ + 1;
#line 260 "sgemlqt.f"
	    ib = min(i__2,i__3);
#line 261 "sgemlqt.f"
	    i__2 = *n - i__ + 1;
#line 261 "sgemlqt.f"
	    slarfb_("R", "T", "F", "R", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 264 "sgemlqt.f"
	}

#line 266 "sgemlqt.f"
    }

#line 268 "sgemlqt.f"
    return 0;

/*     End of SGEMLQT */

} /* sgemlqt_ */

