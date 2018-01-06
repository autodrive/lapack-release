#line 1 "cunm2r.f"
/* cunm2r.f -- translated by f2c (version 20100827).
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

#line 1 "cunm2r.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CUNM2R multiplies a general matrix by the unitary matrix from a QR factorization determined by 
cgeqrf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNM2R + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunm2r.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunm2r.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunm2r.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > CUNM2R overwrites the general complex m-by-n matrix C with */
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
/* > as returned by CGEQRF. Q is of order m if SIDE = 'L' and of order n */
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
/* >          A is COMPLEX array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGEQRF in the first k columns of its array argument A. */
/* >          A is modified by the routine but restored on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by CGEQRF. */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunm2r_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i__, i1, i2, i3, ic, jc, mi, ni, nq;
    static doublecomplex aii;
    static logical left;
    static doublecomplex taui;
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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

#line 200 "cunm2r.f"
    /* Parameter adjustments */
#line 200 "cunm2r.f"
    a_dim1 = *lda;
#line 200 "cunm2r.f"
    a_offset = 1 + a_dim1;
#line 200 "cunm2r.f"
    a -= a_offset;
#line 200 "cunm2r.f"
    --tau;
#line 200 "cunm2r.f"
    c_dim1 = *ldc;
#line 200 "cunm2r.f"
    c_offset = 1 + c_dim1;
#line 200 "cunm2r.f"
    c__ -= c_offset;
#line 200 "cunm2r.f"
    --work;
#line 200 "cunm2r.f"

#line 200 "cunm2r.f"
    /* Function Body */
#line 200 "cunm2r.f"
    *info = 0;
#line 201 "cunm2r.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 202 "cunm2r.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 206 "cunm2r.f"
    if (left) {
#line 207 "cunm2r.f"
	nq = *m;
#line 208 "cunm2r.f"
    } else {
#line 209 "cunm2r.f"
	nq = *n;
#line 210 "cunm2r.f"
    }
#line 211 "cunm2r.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 212 "cunm2r.f"
	*info = -1;
#line 213 "cunm2r.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 214 "cunm2r.f"
	*info = -2;
#line 215 "cunm2r.f"
    } else if (*m < 0) {
#line 216 "cunm2r.f"
	*info = -3;
#line 217 "cunm2r.f"
    } else if (*n < 0) {
#line 218 "cunm2r.f"
	*info = -4;
#line 219 "cunm2r.f"
    } else if (*k < 0 || *k > nq) {
#line 220 "cunm2r.f"
	*info = -5;
#line 221 "cunm2r.f"
    } else if (*lda < max(1,nq)) {
#line 222 "cunm2r.f"
	*info = -7;
#line 223 "cunm2r.f"
    } else if (*ldc < max(1,*m)) {
#line 224 "cunm2r.f"
	*info = -10;
#line 225 "cunm2r.f"
    }
#line 226 "cunm2r.f"
    if (*info != 0) {
#line 227 "cunm2r.f"
	i__1 = -(*info);
#line 227 "cunm2r.f"
	xerbla_("CUNM2R", &i__1, (ftnlen)6);
#line 228 "cunm2r.f"
	return 0;
#line 229 "cunm2r.f"
    }

/*     Quick return if possible */

#line 233 "cunm2r.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 233 "cunm2r.f"
	return 0;
#line 233 "cunm2r.f"
    }

#line 236 "cunm2r.f"
    if (left && ! notran || ! left && notran) {
#line 237 "cunm2r.f"
	i1 = 1;
#line 238 "cunm2r.f"
	i2 = *k;
#line 239 "cunm2r.f"
	i3 = 1;
#line 240 "cunm2r.f"
    } else {
#line 241 "cunm2r.f"
	i1 = *k;
#line 242 "cunm2r.f"
	i2 = 1;
#line 243 "cunm2r.f"
	i3 = -1;
#line 244 "cunm2r.f"
    }

#line 246 "cunm2r.f"
    if (left) {
#line 247 "cunm2r.f"
	ni = *n;
#line 248 "cunm2r.f"
	jc = 1;
#line 249 "cunm2r.f"
    } else {
#line 250 "cunm2r.f"
	mi = *m;
#line 251 "cunm2r.f"
	ic = 1;
#line 252 "cunm2r.f"
    }

#line 254 "cunm2r.f"
    i__1 = i2;
#line 254 "cunm2r.f"
    i__2 = i3;
#line 254 "cunm2r.f"
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 255 "cunm2r.f"
	if (left) {

/*           H(i) or H(i)**H is applied to C(i:m,1:n) */

#line 259 "cunm2r.f"
	    mi = *m - i__ + 1;
#line 260 "cunm2r.f"
	    ic = i__;
#line 261 "cunm2r.f"
	} else {

/*           H(i) or H(i)**H is applied to C(1:m,i:n) */

#line 265 "cunm2r.f"
	    ni = *n - i__ + 1;
#line 266 "cunm2r.f"
	    jc = i__;
#line 267 "cunm2r.f"
	}

/*        Apply H(i) or H(i)**H */

#line 271 "cunm2r.f"
	if (notran) {
#line 272 "cunm2r.f"
	    i__3 = i__;
#line 272 "cunm2r.f"
	    taui.r = tau[i__3].r, taui.i = tau[i__3].i;
#line 273 "cunm2r.f"
	} else {
#line 274 "cunm2r.f"
	    d_cnjg(&z__1, &tau[i__]);
#line 274 "cunm2r.f"
	    taui.r = z__1.r, taui.i = z__1.i;
#line 275 "cunm2r.f"
	}
#line 276 "cunm2r.f"
	i__3 = i__ + i__ * a_dim1;
#line 276 "cunm2r.f"
	aii.r = a[i__3].r, aii.i = a[i__3].i;
#line 277 "cunm2r.f"
	i__3 = i__ + i__ * a_dim1;
#line 277 "cunm2r.f"
	a[i__3].r = 1., a[i__3].i = 0.;
#line 278 "cunm2r.f"
	clarf_(side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &taui, &c__[ic 
		+ jc * c_dim1], ldc, &work[1], (ftnlen)1);
#line 280 "cunm2r.f"
	i__3 = i__ + i__ * a_dim1;
#line 280 "cunm2r.f"
	a[i__3].r = aii.r, a[i__3].i = aii.i;
#line 281 "cunm2r.f"
/* L10: */
#line 281 "cunm2r.f"
    }
#line 282 "cunm2r.f"
    return 0;

/*     End of CUNM2R */

} /* cunm2r_ */

