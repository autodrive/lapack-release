#line 1 "dorm2l.f"
/* dorm2l.f -- translated by f2c (version 20100827).
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

#line 1 "dorm2l.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DORM2L multiplies a general matrix by the orthogonal matrix from a QL factorization determined 
by sgeqlf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORM2L + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2l.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm2l.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm2l.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORM2L overwrites the general real m by n matrix C with */
/* > */
/* >       Q * C  if SIDE = 'L' and TRANS = 'N', or */
/* > */
/* >       Q**T * C  if SIDE = 'L' and TRANS = 'T', or */
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
/* > as returned by DGEQLF. Q is of order m if SIDE = 'L' and of order n */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DGEQLF in the last k columns of its array argument A. */
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
/* >          TAU is DOUBLE PRECISION array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGEQLF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorm2l_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info, ftnlen side_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__, i1, i2, i3, mi, ni, nq;
    static doublereal aii;
    static logical left;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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

#line 200 "dorm2l.f"
    /* Parameter adjustments */
#line 200 "dorm2l.f"
    a_dim1 = *lda;
#line 200 "dorm2l.f"
    a_offset = 1 + a_dim1;
#line 200 "dorm2l.f"
    a -= a_offset;
#line 200 "dorm2l.f"
    --tau;
#line 200 "dorm2l.f"
    c_dim1 = *ldc;
#line 200 "dorm2l.f"
    c_offset = 1 + c_dim1;
#line 200 "dorm2l.f"
    c__ -= c_offset;
#line 200 "dorm2l.f"
    --work;
#line 200 "dorm2l.f"

#line 200 "dorm2l.f"
    /* Function Body */
#line 200 "dorm2l.f"
    *info = 0;
#line 201 "dorm2l.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 202 "dorm2l.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 206 "dorm2l.f"
    if (left) {
#line 207 "dorm2l.f"
	nq = *m;
#line 208 "dorm2l.f"
    } else {
#line 209 "dorm2l.f"
	nq = *n;
#line 210 "dorm2l.f"
    }
#line 211 "dorm2l.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 212 "dorm2l.f"
	*info = -1;
#line 213 "dorm2l.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 214 "dorm2l.f"
	*info = -2;
#line 215 "dorm2l.f"
    } else if (*m < 0) {
#line 216 "dorm2l.f"
	*info = -3;
#line 217 "dorm2l.f"
    } else if (*n < 0) {
#line 218 "dorm2l.f"
	*info = -4;
#line 219 "dorm2l.f"
    } else if (*k < 0 || *k > nq) {
#line 220 "dorm2l.f"
	*info = -5;
#line 221 "dorm2l.f"
    } else if (*lda < max(1,nq)) {
#line 222 "dorm2l.f"
	*info = -7;
#line 223 "dorm2l.f"
    } else if (*ldc < max(1,*m)) {
#line 224 "dorm2l.f"
	*info = -10;
#line 225 "dorm2l.f"
    }
#line 226 "dorm2l.f"
    if (*info != 0) {
#line 227 "dorm2l.f"
	i__1 = -(*info);
#line 227 "dorm2l.f"
	xerbla_("DORM2L", &i__1, (ftnlen)6);
#line 228 "dorm2l.f"
	return 0;
#line 229 "dorm2l.f"
    }

/*     Quick return if possible */

#line 233 "dorm2l.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 233 "dorm2l.f"
	return 0;
#line 233 "dorm2l.f"
    }

#line 236 "dorm2l.f"
    if (left && notran || ! left && ! notran) {
#line 238 "dorm2l.f"
	i1 = 1;
#line 239 "dorm2l.f"
	i2 = *k;
#line 240 "dorm2l.f"
	i3 = 1;
#line 241 "dorm2l.f"
    } else {
#line 242 "dorm2l.f"
	i1 = *k;
#line 243 "dorm2l.f"
	i2 = 1;
#line 244 "dorm2l.f"
	i3 = -1;
#line 245 "dorm2l.f"
    }

#line 247 "dorm2l.f"
    if (left) {
#line 248 "dorm2l.f"
	ni = *n;
#line 249 "dorm2l.f"
    } else {
#line 250 "dorm2l.f"
	mi = *m;
#line 251 "dorm2l.f"
    }

#line 253 "dorm2l.f"
    i__1 = i2;
#line 253 "dorm2l.f"
    i__2 = i3;
#line 253 "dorm2l.f"
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 254 "dorm2l.f"
	if (left) {

/*           H(i) is applied to C(1:m-k+i,1:n) */

#line 258 "dorm2l.f"
	    mi = *m - *k + i__;
#line 259 "dorm2l.f"
	} else {

/*           H(i) is applied to C(1:m,1:n-k+i) */

#line 263 "dorm2l.f"
	    ni = *n - *k + i__;
#line 264 "dorm2l.f"
	}

/*        Apply H(i) */

#line 268 "dorm2l.f"
	aii = a[nq - *k + i__ + i__ * a_dim1];
#line 269 "dorm2l.f"
	a[nq - *k + i__ + i__ * a_dim1] = 1.;
#line 270 "dorm2l.f"
	dlarf_(side, &mi, &ni, &a[i__ * a_dim1 + 1], &c__1, &tau[i__], &c__[
		c_offset], ldc, &work[1], (ftnlen)1);
#line 272 "dorm2l.f"
	a[nq - *k + i__ + i__ * a_dim1] = aii;
#line 273 "dorm2l.f"
/* L10: */
#line 273 "dorm2l.f"
    }
#line 274 "dorm2l.f"
    return 0;

/*     End of DORM2L */

} /* dorm2l_ */

