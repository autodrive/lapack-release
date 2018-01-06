#line 1 "sgemqrt.f"
/* sgemqrt.f -- translated by f2c (version 20100827).
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

#line 1 "sgemqrt.f"
/* > \brief \b SGEMQRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgemqrt
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgemqrt
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgemqrt
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, */
/*                          C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER SIDE, TRANS */
/*       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEMQRT overwrites the general real M-by-N matrix C with */
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
/* > generated using the compact WY representation as returned by SGEQRT. */
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
/* >          = 'T':  Transpose, apply Q**T. */
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
/* >          V is REAL array, dimension (LDV,K) */
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
/* >          T is REAL array, dimension (LDT,K) */
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

/* > \date November 2013 */

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgemqrt_(char *side, char *trans, integer *m, integer *n,
	 integer *k, integer *nb, doublereal *v, integer *ldv, doublereal *t, 
	integer *ldt, doublereal *c__, integer *ldc, doublereal *work, 
	integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, q, ib, kf;
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

#line 205 "sgemqrt.f"
    /* Parameter adjustments */
#line 205 "sgemqrt.f"
    v_dim1 = *ldv;
#line 205 "sgemqrt.f"
    v_offset = 1 + v_dim1;
#line 205 "sgemqrt.f"
    v -= v_offset;
#line 205 "sgemqrt.f"
    t_dim1 = *ldt;
#line 205 "sgemqrt.f"
    t_offset = 1 + t_dim1;
#line 205 "sgemqrt.f"
    t -= t_offset;
#line 205 "sgemqrt.f"
    c_dim1 = *ldc;
#line 205 "sgemqrt.f"
    c_offset = 1 + c_dim1;
#line 205 "sgemqrt.f"
    c__ -= c_offset;
#line 205 "sgemqrt.f"
    --work;
#line 205 "sgemqrt.f"

#line 205 "sgemqrt.f"
    /* Function Body */
#line 205 "sgemqrt.f"
    *info = 0;
#line 206 "sgemqrt.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 207 "sgemqrt.f"
    right = lsame_(side, "R", (ftnlen)1, (ftnlen)1);
#line 208 "sgemqrt.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 209 "sgemqrt.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 211 "sgemqrt.f"
    if (left) {
#line 212 "sgemqrt.f"
	ldwork = max(1,*n);
#line 213 "sgemqrt.f"
	q = *m;
#line 214 "sgemqrt.f"
    } else if (right) {
#line 215 "sgemqrt.f"
	ldwork = max(1,*m);
#line 216 "sgemqrt.f"
	q = *n;
#line 217 "sgemqrt.f"
    }
#line 218 "sgemqrt.f"
    if (! left && ! right) {
#line 219 "sgemqrt.f"
	*info = -1;
#line 220 "sgemqrt.f"
    } else if (! tran && ! notran) {
#line 221 "sgemqrt.f"
	*info = -2;
#line 222 "sgemqrt.f"
    } else if (*m < 0) {
#line 223 "sgemqrt.f"
	*info = -3;
#line 224 "sgemqrt.f"
    } else if (*n < 0) {
#line 225 "sgemqrt.f"
	*info = -4;
#line 226 "sgemqrt.f"
    } else if (*k < 0 || *k > q) {
#line 227 "sgemqrt.f"
	*info = -5;
#line 228 "sgemqrt.f"
    } else if (*nb < 1 || *nb > *k && *k > 0) {
#line 229 "sgemqrt.f"
	*info = -6;
#line 230 "sgemqrt.f"
    } else if (*ldv < max(1,q)) {
#line 231 "sgemqrt.f"
	*info = -8;
#line 232 "sgemqrt.f"
    } else if (*ldt < *nb) {
#line 233 "sgemqrt.f"
	*info = -10;
#line 234 "sgemqrt.f"
    } else if (*ldc < max(1,*m)) {
#line 235 "sgemqrt.f"
	*info = -12;
#line 236 "sgemqrt.f"
    }

#line 238 "sgemqrt.f"
    if (*info != 0) {
#line 239 "sgemqrt.f"
	i__1 = -(*info);
#line 239 "sgemqrt.f"
	xerbla_("SGEMQRT", &i__1, (ftnlen)7);
#line 240 "sgemqrt.f"
	return 0;
#line 241 "sgemqrt.f"
    }

/*     .. Quick return if possible .. */

#line 245 "sgemqrt.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 245 "sgemqrt.f"
	return 0;
#line 245 "sgemqrt.f"
    }

#line 247 "sgemqrt.f"
    if (left && tran) {

#line 249 "sgemqrt.f"
	i__1 = *k;
#line 249 "sgemqrt.f"
	i__2 = *nb;
#line 249 "sgemqrt.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 250 "sgemqrt.f"
	    i__3 = *nb, i__4 = *k - i__ + 1;
#line 250 "sgemqrt.f"
	    ib = min(i__3,i__4);
#line 251 "sgemqrt.f"
	    i__3 = *m - i__ + 1;
#line 251 "sgemqrt.f"
	    slarfb_("L", "T", "F", "C", &i__3, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 254 "sgemqrt.f"
	}

#line 256 "sgemqrt.f"
    } else if (right && notran) {

#line 258 "sgemqrt.f"
	i__2 = *k;
#line 258 "sgemqrt.f"
	i__1 = *nb;
#line 258 "sgemqrt.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 259 "sgemqrt.f"
	    i__3 = *nb, i__4 = *k - i__ + 1;
#line 259 "sgemqrt.f"
	    ib = min(i__3,i__4);
#line 260 "sgemqrt.f"
	    i__3 = *n - i__ + 1;
#line 260 "sgemqrt.f"
	    slarfb_("R", "N", "F", "C", m, &i__3, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 263 "sgemqrt.f"
	}

#line 265 "sgemqrt.f"
    } else if (left && notran) {

#line 267 "sgemqrt.f"
	kf = (*k - 1) / *nb * *nb + 1;
#line 268 "sgemqrt.f"
	i__1 = -(*nb);
#line 268 "sgemqrt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 269 "sgemqrt.f"
	    i__2 = *nb, i__3 = *k - i__ + 1;
#line 269 "sgemqrt.f"
	    ib = min(i__2,i__3);
#line 270 "sgemqrt.f"
	    i__2 = *m - i__ + 1;
#line 270 "sgemqrt.f"
	    slarfb_("L", "N", "F", "C", &i__2, n, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, 
		    &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 273 "sgemqrt.f"
	}

#line 275 "sgemqrt.f"
    } else if (right && tran) {

#line 277 "sgemqrt.f"
	kf = (*k - 1) / *nb * *nb + 1;
#line 278 "sgemqrt.f"
	i__1 = -(*nb);
#line 278 "sgemqrt.f"
	for (i__ = kf; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
#line 279 "sgemqrt.f"
	    i__2 = *nb, i__3 = *k - i__ + 1;
#line 279 "sgemqrt.f"
	    ib = min(i__2,i__3);
#line 280 "sgemqrt.f"
	    i__2 = *n - i__ + 1;
#line 280 "sgemqrt.f"
	    slarfb_("R", "T", "F", "C", m, &i__2, &ib, &v[i__ + i__ * v_dim1],
		     ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 283 "sgemqrt.f"
	}

#line 285 "sgemqrt.f"
    }

#line 287 "sgemqrt.f"
    return 0;

/*     End of SGEMQRT */

} /* sgemqrt_ */

