#line 1 "sorm2l.f"
/* sorm2l.f -- translated by f2c (version 20100827).
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

#line 1 "sorm2l.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SORM2L multiplies a general matrix by the orthogonal matrix from a QL factorization determined 
by sgeqlf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORM2L + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorm2l.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorm2l.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorm2l.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > SORM2L overwrites the general real m by n matrix C with */
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
/* > as returned by SGEQLF. Q is of order m if SIDE = 'L' and of order n */
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
/* >          A is REAL array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          SGEQLF in the last k columns of its array argument A. */
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
/* >          TAU is REAL array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by SGEQLF. */
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
/* Subroutine */ int sorm2l_(char *side, char *trans, integer *m, integer *n, 
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

#line 200 "sorm2l.f"
    /* Parameter adjustments */
#line 200 "sorm2l.f"
    a_dim1 = *lda;
#line 200 "sorm2l.f"
    a_offset = 1 + a_dim1;
#line 200 "sorm2l.f"
    a -= a_offset;
#line 200 "sorm2l.f"
    --tau;
#line 200 "sorm2l.f"
    c_dim1 = *ldc;
#line 200 "sorm2l.f"
    c_offset = 1 + c_dim1;
#line 200 "sorm2l.f"
    c__ -= c_offset;
#line 200 "sorm2l.f"
    --work;
#line 200 "sorm2l.f"

#line 200 "sorm2l.f"
    /* Function Body */
#line 200 "sorm2l.f"
    *info = 0;
#line 201 "sorm2l.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 202 "sorm2l.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     NQ is the order of Q */

#line 206 "sorm2l.f"
    if (left) {
#line 207 "sorm2l.f"
	nq = *m;
#line 208 "sorm2l.f"
    } else {
#line 209 "sorm2l.f"
	nq = *n;
#line 210 "sorm2l.f"
    }
#line 211 "sorm2l.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 212 "sorm2l.f"
	*info = -1;
#line 213 "sorm2l.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 214 "sorm2l.f"
	*info = -2;
#line 215 "sorm2l.f"
    } else if (*m < 0) {
#line 216 "sorm2l.f"
	*info = -3;
#line 217 "sorm2l.f"
    } else if (*n < 0) {
#line 218 "sorm2l.f"
	*info = -4;
#line 219 "sorm2l.f"
    } else if (*k < 0 || *k > nq) {
#line 220 "sorm2l.f"
	*info = -5;
#line 221 "sorm2l.f"
    } else if (*lda < max(1,nq)) {
#line 222 "sorm2l.f"
	*info = -7;
#line 223 "sorm2l.f"
    } else if (*ldc < max(1,*m)) {
#line 224 "sorm2l.f"
	*info = -10;
#line 225 "sorm2l.f"
    }
#line 226 "sorm2l.f"
    if (*info != 0) {
#line 227 "sorm2l.f"
	i__1 = -(*info);
#line 227 "sorm2l.f"
	xerbla_("SORM2L", &i__1, (ftnlen)6);
#line 228 "sorm2l.f"
	return 0;
#line 229 "sorm2l.f"
    }

/*     Quick return if possible */

#line 233 "sorm2l.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 233 "sorm2l.f"
	return 0;
#line 233 "sorm2l.f"
    }

#line 236 "sorm2l.f"
    if (left && notran || ! left && ! notran) {
#line 238 "sorm2l.f"
	i1 = 1;
#line 239 "sorm2l.f"
	i2 = *k;
#line 240 "sorm2l.f"
	i3 = 1;
#line 241 "sorm2l.f"
    } else {
#line 242 "sorm2l.f"
	i1 = *k;
#line 243 "sorm2l.f"
	i2 = 1;
#line 244 "sorm2l.f"
	i3 = -1;
#line 245 "sorm2l.f"
    }

#line 247 "sorm2l.f"
    if (left) {
#line 248 "sorm2l.f"
	ni = *n;
#line 249 "sorm2l.f"
    } else {
#line 250 "sorm2l.f"
	mi = *m;
#line 251 "sorm2l.f"
    }

#line 253 "sorm2l.f"
    i__1 = i2;
#line 253 "sorm2l.f"
    i__2 = i3;
#line 253 "sorm2l.f"
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 254 "sorm2l.f"
	if (left) {

/*           H(i) is applied to C(1:m-k+i,1:n) */

#line 258 "sorm2l.f"
	    mi = *m - *k + i__;
#line 259 "sorm2l.f"
	} else {

/*           H(i) is applied to C(1:m,1:n-k+i) */

#line 263 "sorm2l.f"
	    ni = *n - *k + i__;
#line 264 "sorm2l.f"
	}

/*        Apply H(i) */

#line 268 "sorm2l.f"
	aii = a[nq - *k + i__ + i__ * a_dim1];
#line 269 "sorm2l.f"
	a[nq - *k + i__ + i__ * a_dim1] = 1.;
#line 270 "sorm2l.f"
	slarf_(side, &mi, &ni, &a[i__ * a_dim1 + 1], &c__1, &tau[i__], &c__[
		c_offset], ldc, &work[1], (ftnlen)1);
#line 272 "sorm2l.f"
	a[nq - *k + i__ + i__ * a_dim1] = aii;
#line 273 "sorm2l.f"
/* L10: */
#line 273 "sorm2l.f"
    }
#line 274 "sorm2l.f"
    return 0;

/*     End of SORM2L */

} /* sorm2l_ */

