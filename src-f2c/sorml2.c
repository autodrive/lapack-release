#line 1 "sorml2.f"
/* sorml2.f -- translated by f2c (version 20100827).
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

#line 1 "sorml2.f"
/* > \brief \b SORML2 multiplies a general matrix by the orthogonal matrix from a LQ factorization determined 
by sgelqf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORML2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorml2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorml2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorml2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORML2 overwrites the general real m by n matrix C with */
/* > */
/* >       Q * C  if SIDE = 'L' and TRANS = 'N', or */
/* > */
/* >       Q**T* C  if SIDE = 'L' and TRANS = 'T', or */
/* > */
/* >       C * Q  if SIDE = 'R' and TRANS = 'N', or */
/* > */
/* >       C * Q**T if SIDE = 'R' and TRANS = 'T', */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by SGELQF. Q is of order m if SIDE = 'L' and of order n */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          SGELQF in the first k rows of its array argument A. */
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
/* >          reflector H(i), as returned by SGELQF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
/* >          On entry, the m by n matrix C. */
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

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorml2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info, ftnlen side_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__, i1, i2, i3, ic, jc, mi, ni, nq;
    static doublereal aii;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), xerbla_(char *, integer *, ftnlen);
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

#line 200 "sorml2.f"
    /* Parameter adjustments */
#line 200 "sorml2.f"
    a_dim1 = *lda;
#line 200 "sorml2.f"
    a_offset = 1 + a_dim1;
#line 200 "sorml2.f"
    a -= a_offset;
#line 200 "sorml2.f"
    --tau;
#line 200 "sorml2.f"
    c_dim1 = *ldc;
#line 200 "sorml2.f"
    c_offset = 1 + c_dim1;
#line 200 "sorml2.f"
    c__ -= c_offset;
#line 200 "sorml2.f"
    --work;
#line 200 "sorml2.f"

#line 200 "sorml2.f"
    /* Function Body */
#line 200 "sorml2.f"
    *info = 0;
#line 201 "sorml2.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 202 "sorml2.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 206 "sorml2.f"
    if (left) {
#line 207 "sorml2.f"
	nq = *m;
#line 208 "sorml2.f"
    } else {
#line 209 "sorml2.f"
	nq = *n;
#line 210 "sorml2.f"
    }
#line 211 "sorml2.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 212 "sorml2.f"
	*info = -1;
#line 213 "sorml2.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 214 "sorml2.f"
	*info = -2;
#line 215 "sorml2.f"
    } else if (*m < 0) {
#line 216 "sorml2.f"
	*info = -3;
#line 217 "sorml2.f"
    } else if (*n < 0) {
#line 218 "sorml2.f"
	*info = -4;
#line 219 "sorml2.f"
    } else if (*k < 0 || *k > nq) {
#line 220 "sorml2.f"
	*info = -5;
#line 221 "sorml2.f"
    } else if (*lda < max(1,*k)) {
#line 222 "sorml2.f"
	*info = -7;
#line 223 "sorml2.f"
    } else if (*ldc < max(1,*m)) {
#line 224 "sorml2.f"
	*info = -10;
#line 225 "sorml2.f"
    }
#line 226 "sorml2.f"
    if (*info != 0) {
#line 227 "sorml2.f"
	i__1 = -(*info);
#line 227 "sorml2.f"
	xerbla_("SORML2", &i__1, (ftnlen)6);
#line 228 "sorml2.f"
	return 0;
#line 229 "sorml2.f"
    }

/*     Quick return if possible */

#line 233 "sorml2.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 233 "sorml2.f"
	return 0;
#line 233 "sorml2.f"
    }

#line 236 "sorml2.f"
    if (left && notran || ! left && ! notran) {
#line 238 "sorml2.f"
	i1 = 1;
#line 239 "sorml2.f"
	i2 = *k;
#line 240 "sorml2.f"
	i3 = 1;
#line 241 "sorml2.f"
    } else {
#line 242 "sorml2.f"
	i1 = *k;
#line 243 "sorml2.f"
	i2 = 1;
#line 244 "sorml2.f"
	i3 = -1;
#line 245 "sorml2.f"
    }

#line 247 "sorml2.f"
    if (left) {
#line 248 "sorml2.f"
	ni = *n;
#line 249 "sorml2.f"
	jc = 1;
#line 250 "sorml2.f"
    } else {
#line 251 "sorml2.f"
	mi = *m;
#line 252 "sorml2.f"
	ic = 1;
#line 253 "sorml2.f"
    }

#line 255 "sorml2.f"
    i__1 = i2;
#line 255 "sorml2.f"
    i__2 = i3;
#line 255 "sorml2.f"
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 256 "sorml2.f"
	if (left) {

/*           H(i) is applied to C(i:m,1:n) */

#line 260 "sorml2.f"
	    mi = *m - i__ + 1;
#line 261 "sorml2.f"
	    ic = i__;
#line 262 "sorml2.f"
	} else {

/*           H(i) is applied to C(1:m,i:n) */

#line 266 "sorml2.f"
	    ni = *n - i__ + 1;
#line 267 "sorml2.f"
	    jc = i__;
#line 268 "sorml2.f"
	}

/*        Apply H(i) */

#line 272 "sorml2.f"
	aii = a[i__ + i__ * a_dim1];
#line 273 "sorml2.f"
	a[i__ + i__ * a_dim1] = 1.;
#line 274 "sorml2.f"
	slarf_(side, &mi, &ni, &a[i__ + i__ * a_dim1], lda, &tau[i__], &c__[
		ic + jc * c_dim1], ldc, &work[1], (ftnlen)1);
#line 276 "sorml2.f"
	a[i__ + i__ * a_dim1] = aii;
#line 277 "sorml2.f"
/* L10: */
#line 277 "sorml2.f"
    }
#line 278 "sorml2.f"
    return 0;

/*     End of SORML2 */

} /* sorml2_ */

