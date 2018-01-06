#line 1 "zunmr3.f"
/* zunmr3.f -- translated by f2c (version 20100827).
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

#line 1 "zunmr3.f"
/* > \brief \b ZUNMR3 multiplies a general matrix by the unitary matrix from a RZ factorization determined by 
ctzrzf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNMR3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmr3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmr3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmr3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, L, LDA, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNMR3 overwrites the general complex m by n matrix C with */
/* > */
/* >       Q * C  if SIDE = 'L' and TRANS = 'N', or */
/* > */
/* >       Q**H* C  if SIDE = 'L' and TRANS = 'C', or */
/* > */
/* >       C * Q  if SIDE = 'R' and TRANS = 'N', or */
/* > */
/* >       C * Q**H if SIDE = 'R' and TRANS = 'C', */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by ZTZRZF. Q is of order m if SIDE = 'L' and of order n */
/* > if SIDE = 'R'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left */
/* >          = 'R': apply Q or Q**H from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply Q  (No transpose) */
/* >          = 'C': apply Q**H (Conjugate transpose) */
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
/* >          A is COMPLEX*16 array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          ZTZRZF in the last k rows of its array argument A. */
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
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZTZRZF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the m-by-n matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
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
/* >          WORK is COMPLEX*16 array, dimension */
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

/* > \ingroup complex16OTHERcomputational */

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
/* Subroutine */ int zunmr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, i1, i2, i3, ja, ic, jc, mi, ni, nq;
    static logical left;
    static doublecomplex taui;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarz_(char *, integer *, integer *, integer *
	    , doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
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

#line 215 "zunmr3.f"
    /* Parameter adjustments */
#line 215 "zunmr3.f"
    a_dim1 = *lda;
#line 215 "zunmr3.f"
    a_offset = 1 + a_dim1;
#line 215 "zunmr3.f"
    a -= a_offset;
#line 215 "zunmr3.f"
    --tau;
#line 215 "zunmr3.f"
    c_dim1 = *ldc;
#line 215 "zunmr3.f"
    c_offset = 1 + c_dim1;
#line 215 "zunmr3.f"
    c__ -= c_offset;
#line 215 "zunmr3.f"
    --work;
#line 215 "zunmr3.f"

#line 215 "zunmr3.f"
    /* Function Body */
#line 215 "zunmr3.f"
    *info = 0;
#line 216 "zunmr3.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 217 "zunmr3.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 221 "zunmr3.f"
    if (left) {
#line 222 "zunmr3.f"
	nq = *m;
#line 223 "zunmr3.f"
    } else {
#line 224 "zunmr3.f"
	nq = *n;
#line 225 "zunmr3.f"
    }
#line 226 "zunmr3.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 227 "zunmr3.f"
	*info = -1;
#line 228 "zunmr3.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 229 "zunmr3.f"
	*info = -2;
#line 230 "zunmr3.f"
    } else if (*m < 0) {
#line 231 "zunmr3.f"
	*info = -3;
#line 232 "zunmr3.f"
    } else if (*n < 0) {
#line 233 "zunmr3.f"
	*info = -4;
#line 234 "zunmr3.f"
    } else if (*k < 0 || *k > nq) {
#line 235 "zunmr3.f"
	*info = -5;
#line 236 "zunmr3.f"
    } else if (*l < 0 || left && *l > *m || ! left && *l > *n) {
#line 238 "zunmr3.f"
	*info = -6;
#line 239 "zunmr3.f"
    } else if (*lda < max(1,*k)) {
#line 240 "zunmr3.f"
	*info = -8;
#line 241 "zunmr3.f"
    } else if (*ldc < max(1,*m)) {
#line 242 "zunmr3.f"
	*info = -11;
#line 243 "zunmr3.f"
    }
#line 244 "zunmr3.f"
    if (*info != 0) {
#line 245 "zunmr3.f"
	i__1 = -(*info);
#line 245 "zunmr3.f"
	xerbla_("ZUNMR3", &i__1, (ftnlen)6);
#line 246 "zunmr3.f"
	return 0;
#line 247 "zunmr3.f"
    }

/*     Quick return if possible */

#line 251 "zunmr3.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 251 "zunmr3.f"
	return 0;
#line 251 "zunmr3.f"
    }

#line 254 "zunmr3.f"
    if (left && ! notran || ! left && notran) {
#line 255 "zunmr3.f"
	i1 = 1;
#line 256 "zunmr3.f"
	i2 = *k;
#line 257 "zunmr3.f"
	i3 = 1;
#line 258 "zunmr3.f"
    } else {
#line 259 "zunmr3.f"
	i1 = *k;
#line 260 "zunmr3.f"
	i2 = 1;
#line 261 "zunmr3.f"
	i3 = -1;
#line 262 "zunmr3.f"
    }

#line 264 "zunmr3.f"
    if (left) {
#line 265 "zunmr3.f"
	ni = *n;
#line 266 "zunmr3.f"
	ja = *m - *l + 1;
#line 267 "zunmr3.f"
	jc = 1;
#line 268 "zunmr3.f"
    } else {
#line 269 "zunmr3.f"
	mi = *m;
#line 270 "zunmr3.f"
	ja = *n - *l + 1;
#line 271 "zunmr3.f"
	ic = 1;
#line 272 "zunmr3.f"
    }

#line 274 "zunmr3.f"
    i__1 = i2;
#line 274 "zunmr3.f"
    i__2 = i3;
#line 274 "zunmr3.f"
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 275 "zunmr3.f"
	if (left) {

/*           H(i) or H(i)**H is applied to C(i:m,1:n) */

#line 279 "zunmr3.f"
	    mi = *m - i__ + 1;
#line 280 "zunmr3.f"
	    ic = i__;
#line 281 "zunmr3.f"
	} else {

/*           H(i) or H(i)**H is applied to C(1:m,i:n) */

#line 285 "zunmr3.f"
	    ni = *n - i__ + 1;
#line 286 "zunmr3.f"
	    jc = i__;
#line 287 "zunmr3.f"
	}

/*        Apply H(i) or H(i)**H */

#line 291 "zunmr3.f"
	if (notran) {
#line 292 "zunmr3.f"
	    i__3 = i__;
#line 292 "zunmr3.f"
	    taui.r = tau[i__3].r, taui.i = tau[i__3].i;
#line 293 "zunmr3.f"
	} else {
#line 294 "zunmr3.f"
	    d_cnjg(&z__1, &tau[i__]);
#line 294 "zunmr3.f"
	    taui.r = z__1.r, taui.i = z__1.i;
#line 295 "zunmr3.f"
	}
#line 296 "zunmr3.f"
	zlarz_(side, &mi, &ni, l, &a[i__ + ja * a_dim1], lda, &taui, &c__[ic 
		+ jc * c_dim1], ldc, &work[1], (ftnlen)1);

#line 299 "zunmr3.f"
/* L10: */
#line 299 "zunmr3.f"
    }

#line 301 "zunmr3.f"
    return 0;

/*     End of ZUNMR3 */

} /* zunmr3_ */

