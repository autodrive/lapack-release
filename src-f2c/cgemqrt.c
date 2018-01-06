#line 1 "cgemqrt.f"
/* cgemqrt.f -- translated by f2c (version 20100827).
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

#line 1 "cgemqrt.f"
/* > \brief \b CGEMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgemqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgemqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgemqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, */
/*                          C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEMQRT overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q C            C Q */
/* > TRANS = 'C':    Q**H C            C Q**H */
/* > */
/* > where Q is a complex orthogonal matrix defined as the product of K */
/* > elementary reflectors: */
/* > */
/* >       Q = H(1) H(2) . . . H(K) = I - V T V**H */
/* > */
/* > generated using the compact WY representation as returned by CGEQRT. */
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
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The block size used for the storage of T.  K >= NB >= 1. */
/* >          This must be the same value of NB used to generate T */
/* >          in CGEQRT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (LDV,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGEQRT in the first K columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,K) */
/* >          The upper triangular factors of the block reflectors */
/* >          as returned by CGEQRT, stored as a NB-by-N matrix. */
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
/* >          C is COMPLEX array, dimension (LDC,N) */
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
/* >          WORK is COMPLEX array. The dimension of WORK is */
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

/* > \date November 2013 */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgemqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *nb, doublecomplex *v, integer *ldv, 
	doublecomplex *t, integer *ldt, doublecomplex *c__, integer *ldc, 
	doublecomplex *work, integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, q, ib, kf;
    static logical left, tran;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int clarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static integer ldwork;


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

#line 205 "cgemqrt.f"
    /* Parameter adjustments */
#line 205 "cgemqrt.f"
    v_dim1 = *ldv;
#line 205 "cgemqrt.f"
    v_offset = 1 + v_dim1;
#line 205 "cgemqrt.f"
    v -= v_offset;
#line 205 "cgemqrt.f"
    t_dim1 = *ldt;
#line 205 "cgemqrt.f"
    t_offset = 1 + t_dim1;
#line 205 "cgemqrt.f"
    t -= t_offset;
#line 205 "cgemqrt.f"
    c_dim1 = *ldc;
#line 205 "cgemqrt.f"
    c_offset = 1 + c_dim1;
#line 205 "cgemqrt.f"
    c__ -= c_offset;
#line 205 "cgemqrt.f"
    --work;
#line 205 "cgemqrt.f"

#line 205 "cgemqrt.f"
    /* Function Body */
#line 205 "cgemqrt.f"
    *info = 0;
#line 206 "cgemqrt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 207 "cgemqrt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 208 "cgemqrt.f"
    tran = lsame_(trans, "C", (ftnlen)1, (ftnlen)1);
#line 209 "cgemqrt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 211 "cgemqrt.f"
    if (left) {
#line 212 "cgemqrt.f"
	ldwork = max(1,*n);
#line 213 "cgemqrt.f"
	q = *m;
#line 214 "cgemqrt.f"
    } else if (right) {
#line 215 "cgemqrt.f"
	ldwork = max(1,*m);
#line 216 "cgemqrt.f"
	q = *n;
#line 217 "cgemqrt.f"
    }
#line 218 "cgemqrt.f"
    if (! left && ! right) {
#line 219 "cgemqrt.f"
	*info = -1;
#line 220 "cgemqrt.f"
    } else if (! tran && ! notran) {
#line 221 "cgemqrt.f"
	*info = -2;
#line 222 "cgemqrt.f"
    } else if (*m < 0) {
#line 223 "cgemqrt.f"
	*info = -3;
#line 224 "cgemqrt.f"
    } else if (*n < 0) {
#line 225 "cgemqrt.f"
	*info = -4;
#line 226 "cgemqrt.f"
    } else if (*k < 0 || *k > q) {
#line 227 "cgemqrt.f"
	*info = -5;
#line 228 "cgemqrt.f"
    } else if (*nb < 1 || *nb > *k && *k > 0) {
#line 229 "cgemqrt.f"
	*info = -6;
#line 230 "cgemqrt.f"
    } else if (*ldv < max(1,q)) {
#line 231 "cgemqrt.f"
	*info = -8;
#line 232 "cgemqrt.f"
    } else if (*ldt < *nb) {
#line 233 "cgemqrt.f"
	*info = -10;
#line 234 "cgemqrt.f"
    } else if (*ldc < max(1,*m)) {
#line 235 "cgemqrt.f"
	*info = -12;
#line 236 "cgemqrt.f"
    }

#line 238 "cgemqrt.f"
    if (*info != 0) {
#line 239 "cgemqrt.f"
	i__1 = -(*info);
#line 239 "cgemqrt.f"
	xerbla_("CGEMQRT", &i__1, (ftnlen)7);
#line 240 "cgemqrt.f"
	return 0;
#line 241 "cgemqrt.f"
    }

/*     .. Quick return if possible .. */

#line 245 "cgemqrt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 245 "cgemqrt.f"
	return 0;
#line 245 "cgemqrt.f"
    }

#line 247 "cgemqrt.f"
    if (left && tran) {

#line 249 "cgemqrt.f"
	i__1 = *k;
#line 249 "cgemqrt.f"
	i__2 = *nb;
#line 249 "cgemqrt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 250 "cgemqrt.f"
	    i__3 = *nb, i__4 = *k - i__ + 1;
#line 250 "cgemqrt.f"
	    ib = min(i__3,i__4);
#line 251 "cgemqrt.f"
	    i__3 = *m - i__ + 1;
#line 251 "cgemqrt.f"
	    clarfb_("L", "C", "F", "C", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 254 "cgemqrt.f"
	}

#line 256 "cgemqrt.f"
    } else if (right && notran) {

#line 258 "cgemqrt.f"
	i__2 = *k;
#line 258 "cgemqrt.f"
	i__1 = *nb;
#line 258 "cgemqrt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 259 "cgemqrt.f"
	    i__3 = *nb, i__4 = *k - i__ + 1;
#line 259 "cgemqrt.f"
	    ib = min(i__3,i__4);
#line 260 "cgemqrt.f"
	    i__3 = *n - i__ + 1;
#line 260 "cgemqrt.f"
	    clarfb_("R", "N", "F", "C", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 263 "cgemqrt.f"
	}

#line 265 "cgemqrt.f"
    } else if (left && notran) {

#line 267 "cgemqrt.f"
	kf = (*k - 1) / *nb * *nb + 1;
#line 268 "cgemqrt.f"
	i__1 = -(*nb);
#line 268 "cgemqrt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 269 "cgemqrt.f"
	    i__2 = *nb, i__3 = *k - i__ + 1;
#line 269 "cgemqrt.f"
	    ib = min(i__2,i__3);
#line 270 "cgemqrt.f"
	    i__2 = *m - i__ + 1;
#line 270 "cgemqrt.f"
	    clarfb_("L", "N", "F", "C", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 273 "cgemqrt.f"
	}

#line 275 "cgemqrt.f"
    } else if (right && tran) {

#line 277 "cgemqrt.f"
	kf = (*k - 1) / *nb * *nb + 1;
#line 278 "cgemqrt.f"
	i__1 = -(*nb);
#line 278 "cgemqrt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 279 "cgemqrt.f"
	    i__2 = *nb, i__3 = *k - i__ + 1;
#line 279 "cgemqrt.f"
	    ib = min(i__2,i__3);
#line 280 "cgemqrt.f"
	    i__2 = *n - i__ + 1;
#line 280 "cgemqrt.f"
	    clarfb_("R", "C", "F", "C", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 283 "cgemqrt.f"
	}

#line 285 "cgemqrt.f"
    }

#line 287 "cgemqrt.f"
    return 0;

/*     End of CGEMQRT */

} /* cgemqrt_ */

