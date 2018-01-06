#line 1 "sormr3.f"
/* sormr3.f -- translated by f2c (version 20100827).
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

#line 1 "sormr3.f"
/* > \brief \b SORMR3 multiplies a general matrix by the orthogonal matrix from a RZ factorization determined 
by stzrzf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORMR3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormr3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormr3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormr3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, L, LDA, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORMR3 overwrites the general real m by n matrix C with */
/* > */
/* >       Q * C  if SIDE = 'L' and TRANS = 'N', or */
/* > */
/* >       Q**T* C  if SIDE = 'L' and TRANS = 'C', or */
/* > */
/* >       C * Q  if SIDE = 'R' and TRANS = 'N', or */
/* > */
/* >       C * Q**T if SIDE = 'R' and TRANS = 'C', */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by STZRZF. Q is of order m if SIDE = 'L' and of order n */
/* > if SIDE = 'R'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**T from the Left */
/* >          = 'R': apply Q or Q**T from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply Q  (No transpose) */
/* >          = 'T': apply Q**T (Transpose) */
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
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of columns of the matrix A containing */
/* >          the meaningful part of the Householder reflectors. */
/* >          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          STZRZF in the last k rows of its array argument A. */
/* >          A is modified by the routine but restored on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,K). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by STZRZF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
/* >          On entry, the m-by-n matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
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
/* >          WORK is REAL array, dimension */
/* >                                   (N) if SIDE = 'L', */
/* >                                   (M) if SIDE = 'R' */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sormr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *info, 
	ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__, i1, i2, i3, ja, ic, jc, mi, ni, nq;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarz_(char *, integer *, integer *, integer *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical notran;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 214 "sormr3.f"
    /* Parameter adjustments */
#line 214 "sormr3.f"
    a_dim1 = *lda;
#line 214 "sormr3.f"
    a_offset = 1 + a_dim1;
#line 214 "sormr3.f"
    a -= a_offset;
#line 214 "sormr3.f"
    --tau;
#line 214 "sormr3.f"
    c_dim1 = *ldc;
#line 214 "sormr3.f"
    c_offset = 1 + c_dim1;
#line 214 "sormr3.f"
    c__ -= c_offset;
#line 214 "sormr3.f"
    --work;
#line 214 "sormr3.f"

#line 214 "sormr3.f"
    /* Function Body */
#line 214 "sormr3.f"
    *info = 0;
#line 215 "sormr3.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 216 "sormr3.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 220 "sormr3.f"
    if (left) {
#line 221 "sormr3.f"
	nq = *m;
#line 222 "sormr3.f"
    } else {
#line 223 "sormr3.f"
	nq = *n;
#line 224 "sormr3.f"
    }
#line 225 "sormr3.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 226 "sormr3.f"
	*info = -1;
#line 227 "sormr3.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 228 "sormr3.f"
	*info = -2;
#line 229 "sormr3.f"
    } else if (*m < 0) {
#line 230 "sormr3.f"
	*info = -3;
#line 231 "sormr3.f"
    } else if (*n < 0) {
#line 232 "sormr3.f"
	*info = -4;
#line 233 "sormr3.f"
    } else if (*k < 0 || *k > nq) {
#line 234 "sormr3.f"
	*info = -5;
#line 235 "sormr3.f"
    } else if (*l < 0 || left && *l > *m || ! left && *l > *n) {
#line 237 "sormr3.f"
	*info = -6;
#line 238 "sormr3.f"
    } else if (*lda < max(1,*k)) {
#line 239 "sormr3.f"
	*info = -8;
#line 240 "sormr3.f"
    } else if (*ldc < max(1,*m)) {
#line 241 "sormr3.f"
	*info = -11;
#line 242 "sormr3.f"
    }
#line 243 "sormr3.f"
    if (*info != 0) {
#line 244 "sormr3.f"
	i__1 = -(*info);
#line 244 "sormr3.f"
	xerbla_("SORMR3", &i__1, (ftnlen)6);
#line 245 "sormr3.f"
	return 0;
#line 246 "sormr3.f"
    }

/*     Quick return if possible */

#line 250 "sormr3.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 250 "sormr3.f"
	return 0;
#line 250 "sormr3.f"
    }

#line 253 "sormr3.f"
    if (left && ! notran || ! left && notran) {
#line 254 "sormr3.f"
	i1 = 1;
#line 255 "sormr3.f"
	i2 = *k;
#line 256 "sormr3.f"
	i3 = 1;
#line 257 "sormr3.f"
    } else {
#line 258 "sormr3.f"
	i1 = *k;
#line 259 "sormr3.f"
	i2 = 1;
#line 260 "sormr3.f"
	i3 = -1;
#line 261 "sormr3.f"
    }

#line 263 "sormr3.f"
    if (left) {
#line 264 "sormr3.f"
	ni = *n;
#line 265 "sormr3.f"
	ja = *m - *l + 1;
#line 266 "sormr3.f"
	jc = 1;
#line 267 "sormr3.f"
    } else {
#line 268 "sormr3.f"
	mi = *m;
#line 269 "sormr3.f"
	ja = *n - *l + 1;
#line 270 "sormr3.f"
	ic = 1;
#line 271 "sormr3.f"
    }

#line 273 "sormr3.f"
    i__1 = i2;
#line 273 "sormr3.f"
    i__2 = i3;
#line 273 "sormr3.f"
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 274 "sormr3.f"
	if (left) {

/*           H(i) or H(i)**T is applied to C(i:m,1:n) */

#line 278 "sormr3.f"
	    mi = *m - i__ + 1;
#line 279 "sormr3.f"
	    ic = i__;
#line 280 "sormr3.f"
	} else {

/*           H(i) or H(i)**T is applied to C(1:m,i:n) */

#line 284 "sormr3.f"
	    ni = *n - i__ + 1;
#line 285 "sormr3.f"
	    jc = i__;
#line 286 "sormr3.f"
	}

/*        Apply H(i) or H(i)**T */

#line 290 "sormr3.f"
	slarz_(side, &mi, &ni, l, &a[i__ + ja * a_dim1], lda, &tau[i__], &c__[
		ic + jc * c_dim1], ldc, &work[1], (ftnlen)1);

#line 293 "sormr3.f"
/* L10: */
#line 293 "sormr3.f"
    }

#line 295 "sormr3.f"
    return 0;

/*     End of SORMR3 */

} /* sormr3_ */

