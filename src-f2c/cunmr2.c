#line 1 "cunmr2.f"
/* cunmr2.f -- translated by f2c (version 20100827).
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

#line 1 "cunmr2.f"
/* > \brief \b CUNMR2 multiplies a general matrix by the unitary matrix from a RQ factorization determined by 
cgerqf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMR2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmr2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmr2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmr2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNMR2 overwrites the general complex m-by-n matrix C with */
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
/* >       Q = H(1)**H H(2)**H . . . H(k)**H */
/* > */
/* > as returned by CGERQF. Q is of order m if SIDE = 'L' and of order n */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGERQF in the last k rows of its array argument A. */
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
/* >          TAU is COMPLEX array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by CGERQF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC,N) */
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
/* >          WORK is COMPLEX array, dimension */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunmr2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info, 
	ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, i1, i2, i3, mi, ni, nq;
    static doublecomplex aii;
    static logical left;
    static doublecomplex taui;
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
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

/*     .. Parameters .. */
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

/*     Test the input arguments */

#line 200 "cunmr2.f"
    /* Parameter adjustments */
#line 200 "cunmr2.f"
    a_dim1 = *lda;
#line 200 "cunmr2.f"
    a_offset = 1 + a_dim1;
#line 200 "cunmr2.f"
    a -= a_offset;
#line 200 "cunmr2.f"
    --tau;
#line 200 "cunmr2.f"
    c_dim1 = *ldc;
#line 200 "cunmr2.f"
    c_offset = 1 + c_dim1;
#line 200 "cunmr2.f"
    c__ -= c_offset;
#line 200 "cunmr2.f"
    --work;
#line 200 "cunmr2.f"

#line 200 "cunmr2.f"
    /* Function Body */
#line 200 "cunmr2.f"
    *info = 0;
#line 201 "cunmr2.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 202 "cunmr2.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 206 "cunmr2.f"
    if (left) {
#line 207 "cunmr2.f"
	nq = *m;
#line 208 "cunmr2.f"
    } else {
#line 209 "cunmr2.f"
	nq = *n;
#line 210 "cunmr2.f"
    }
#line 211 "cunmr2.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 212 "cunmr2.f"
	*info = -1;
#line 213 "cunmr2.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 214 "cunmr2.f"
	*info = -2;
#line 215 "cunmr2.f"
    } else if (*m < 0) {
#line 216 "cunmr2.f"
	*info = -3;
#line 217 "cunmr2.f"
    } else if (*n < 0) {
#line 218 "cunmr2.f"
	*info = -4;
#line 219 "cunmr2.f"
    } else if (*k < 0 || *k > nq) {
#line 220 "cunmr2.f"
	*info = -5;
#line 221 "cunmr2.f"
    } else if (*lda < max(1,*k)) {
#line 222 "cunmr2.f"
	*info = -7;
#line 223 "cunmr2.f"
    } else if (*ldc < max(1,*m)) {
#line 224 "cunmr2.f"
	*info = -10;
#line 225 "cunmr2.f"
    }
#line 226 "cunmr2.f"
    if (*info != 0) {
#line 227 "cunmr2.f"
	i__1 = -(*info);
#line 227 "cunmr2.f"
	xerbla_("CUNMR2", &i__1, (ftnlen)6);
#line 228 "cunmr2.f"
	return 0;
#line 229 "cunmr2.f"
    }

/*     Quick return if possible */

#line 233 "cunmr2.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 233 "cunmr2.f"
	return 0;
#line 233 "cunmr2.f"
    }

#line 236 "cunmr2.f"
    if (left && ! notran || ! left && notran) {
#line 237 "cunmr2.f"
	i1 = 1;
#line 238 "cunmr2.f"
	i2 = *k;
#line 239 "cunmr2.f"
	i3 = 1;
#line 240 "cunmr2.f"
    } else {
#line 241 "cunmr2.f"
	i1 = *k;
#line 242 "cunmr2.f"
	i2 = 1;
#line 243 "cunmr2.f"
	i3 = -1;
#line 244 "cunmr2.f"
    }

#line 246 "cunmr2.f"
    if (left) {
#line 247 "cunmr2.f"
	ni = *n;
#line 248 "cunmr2.f"
    } else {
#line 249 "cunmr2.f"
	mi = *m;
#line 250 "cunmr2.f"
    }

#line 252 "cunmr2.f"
    i__1 = i2;
#line 252 "cunmr2.f"
    i__2 = i3;
#line 252 "cunmr2.f"
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 253 "cunmr2.f"
	if (left) {

/*           H(i) or H(i)**H is applied to C(1:m-k+i,1:n) */

#line 257 "cunmr2.f"
	    mi = *m - *k + i__;
#line 258 "cunmr2.f"
	} else {

/*           H(i) or H(i)**H is applied to C(1:m,1:n-k+i) */

#line 262 "cunmr2.f"
	    ni = *n - *k + i__;
#line 263 "cunmr2.f"
	}

/*        Apply H(i) or H(i)**H */

#line 267 "cunmr2.f"
	if (notran) {
#line 268 "cunmr2.f"
	    d_cnjg(&z__1, &tau[i__]);
#line 268 "cunmr2.f"
	    taui.r = z__1.r, taui.i = z__1.i;
#line 269 "cunmr2.f"
	} else {
#line 270 "cunmr2.f"
	    i__3 = i__;
#line 270 "cunmr2.f"
	    taui.r = tau[i__3].r, taui.i = tau[i__3].i;
#line 271 "cunmr2.f"
	}
#line 272 "cunmr2.f"
	i__3 = nq - *k + i__ - 1;
#line 272 "cunmr2.f"
	clacgv_(&i__3, &a[i__ + a_dim1], lda);
#line 273 "cunmr2.f"
	i__3 = i__ + (nq - *k + i__) * a_dim1;
#line 273 "cunmr2.f"
	aii.r = a[i__3].r, aii.i = a[i__3].i;
#line 274 "cunmr2.f"
	i__3 = i__ + (nq - *k + i__) * a_dim1;
#line 274 "cunmr2.f"
	a[i__3].r = 1., a[i__3].i = 0.;
#line 275 "cunmr2.f"
	clarf_(side, &mi, &ni, &a[i__ + a_dim1], lda, &taui, &c__[c_offset], 
		ldc, &work[1], (ftnlen)1);
#line 276 "cunmr2.f"
	i__3 = i__ + (nq - *k + i__) * a_dim1;
#line 276 "cunmr2.f"
	a[i__3].r = aii.r, a[i__3].i = aii.i;
#line 277 "cunmr2.f"
	i__3 = nq - *k + i__ - 1;
#line 277 "cunmr2.f"
	clacgv_(&i__3, &a[i__ + a_dim1], lda);
#line 278 "cunmr2.f"
/* L10: */
#line 278 "cunmr2.f"
    }
#line 279 "cunmr2.f"
    return 0;

/*     End of CUNMR2 */

} /* cunmr2_ */

