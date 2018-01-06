#line 1 "dormrq.f"
/* dormrq.f -- translated by f2c (version 20100827).
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

#line 1 "dormrq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b DORMRQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORMRQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormrq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormrq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormrq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORMRQ overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by DGERQF. Q is of order M if SIDE = 'L' and of order N */
/* > if SIDE = 'R'. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DGERQF in the last k rows of its array argument A. */
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
/* >          TAU is DOUBLE PRECISION array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGERQF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If SIDE = 'L', LWORK >= max(1,N); */
/* >          if SIDE = 'R', LWORK >= max(1,M). */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \date November 2015 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dormrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info, 
	ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2], i__4, 
	    i__5;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, i1, i2, i3, ib, nb, mi, ni, nq, nw, iwt;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int dormr2_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen), dlarfb_(char 
	    *, char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), dlarft_(char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical notran;
    static integer ldwork;
    static char transt[1];
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 211 "dormrq.f"
    /* Parameter adjustments */
#line 211 "dormrq.f"
    a_dim1 = *lda;
#line 211 "dormrq.f"
    a_offset = 1 + a_dim1;
#line 211 "dormrq.f"
    a -= a_offset;
#line 211 "dormrq.f"
    --tau;
#line 211 "dormrq.f"
    c_dim1 = *ldc;
#line 211 "dormrq.f"
    c_offset = 1 + c_dim1;
#line 211 "dormrq.f"
    c__ -= c_offset;
#line 211 "dormrq.f"
    --work;
#line 211 "dormrq.f"

#line 211 "dormrq.f"
    /* Function Body */
#line 211 "dormrq.f"
    *info = 0;
#line 212 "dormrq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 213 "dormrq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 214 "dormrq.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 218 "dormrq.f"
    if (left) {
#line 219 "dormrq.f"
	nq = *m;
#line 220 "dormrq.f"
	nw = max(1,*n);
#line 221 "dormrq.f"
    } else {
#line 222 "dormrq.f"
	nq = *n;
#line 223 "dormrq.f"
	nw = max(1,*m);
#line 224 "dormrq.f"
    }
#line 225 "dormrq.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 226 "dormrq.f"
	*info = -1;
#line 227 "dormrq.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 228 "dormrq.f"
	*info = -2;
#line 229 "dormrq.f"
    } else if (*m < 0) {
#line 230 "dormrq.f"
	*info = -3;
#line 231 "dormrq.f"
    } else if (*n < 0) {
#line 232 "dormrq.f"
	*info = -4;
#line 233 "dormrq.f"
    } else if (*k < 0 || *k > nq) {
#line 234 "dormrq.f"
	*info = -5;
#line 235 "dormrq.f"
    } else if (*lda < max(1,*k)) {
#line 236 "dormrq.f"
	*info = -7;
#line 237 "dormrq.f"
    } else if (*ldc < max(1,*m)) {
#line 238 "dormrq.f"
	*info = -10;
#line 239 "dormrq.f"
    } else if (*lwork < nw && ! lquery) {
#line 240 "dormrq.f"
	*info = -12;
#line 241 "dormrq.f"
    }

#line 243 "dormrq.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 247 "dormrq.f"
	if (*m == 0 || *n == 0) {
#line 248 "dormrq.f"
	    lwkopt = 1;
#line 249 "dormrq.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 250 "dormrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 250 "dormrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 250 "dormrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 250 "dormrq.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 250 "dormrq.f"
	    nb = min(i__1,i__2);
#line 252 "dormrq.f"
	    lwkopt = nw * nb + 4160;
#line 253 "dormrq.f"
	}
#line 254 "dormrq.f"
	work[1] = (doublereal) lwkopt;
#line 255 "dormrq.f"
    }

#line 257 "dormrq.f"
    if (*info != 0) {
#line 258 "dormrq.f"
	i__1 = -(*info);
#line 258 "dormrq.f"
	xerbla_("DORMRQ", &i__1, (ftnlen)6);
#line 259 "dormrq.f"
	return 0;
#line 260 "dormrq.f"
    } else if (lquery) {
#line 261 "dormrq.f"
	return 0;
#line 262 "dormrq.f"
    }

/*     Quick return if possible */

#line 266 "dormrq.f"
    if (*m == 0 || *n == 0) {
#line 267 "dormrq.f"
	return 0;
#line 268 "dormrq.f"
    }

#line 270 "dormrq.f"
    nbmin = 2;
#line 271 "dormrq.f"
    ldwork = nw;
#line 272 "dormrq.f"
    if (nb > 1 && nb < *k) {
#line 273 "dormrq.f"
	if (*lwork < nw * nb + 4160) {
#line 274 "dormrq.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 275 "dormrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 275 "dormrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 275 "dormrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 275 "dormrq.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 275 "dormrq.f"
	    nbmin = max(i__1,i__2);
#line 277 "dormrq.f"
	}
#line 278 "dormrq.f"
    }

#line 280 "dormrq.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 284 "dormrq.f"
	dormr2_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 286 "dormrq.f"
    } else {

/*        Use blocked code */

#line 290 "dormrq.f"
	iwt = nw * nb + 1;
#line 291 "dormrq.f"
	if (left && ! notran || ! left && notran) {
#line 293 "dormrq.f"
	    i1 = 1;
#line 294 "dormrq.f"
	    i2 = *k;
#line 295 "dormrq.f"
	    i3 = nb;
#line 296 "dormrq.f"
	} else {
#line 297 "dormrq.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 298 "dormrq.f"
	    i2 = 1;
#line 299 "dormrq.f"
	    i3 = -nb;
#line 300 "dormrq.f"
	}

#line 302 "dormrq.f"
	if (left) {
#line 303 "dormrq.f"
	    ni = *n;
#line 304 "dormrq.f"
	} else {
#line 305 "dormrq.f"
	    mi = *m;
#line 306 "dormrq.f"
	}

#line 308 "dormrq.f"
	if (notran) {
#line 309 "dormrq.f"
	    *(unsigned char *)transt = 'T';
#line 310 "dormrq.f"
	} else {
#line 311 "dormrq.f"
	    *(unsigned char *)transt = 'N';
#line 312 "dormrq.f"
	}

#line 314 "dormrq.f"
	i__1 = i2;
#line 314 "dormrq.f"
	i__2 = i3;
#line 314 "dormrq.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 315 "dormrq.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 315 "dormrq.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 320 "dormrq.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 320 "dormrq.f"
	    dlarft_("Backward", "Rowwise", &i__4, &ib, &a[i__ + a_dim1], lda, 
		    &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)7);
#line 322 "dormrq.f"
	    if (left) {

/*              H or H**T is applied to C(1:m-k+i+ib-1,1:n) */

#line 326 "dormrq.f"
		mi = *m - *k + i__ + ib - 1;
#line 327 "dormrq.f"
	    } else {

/*              H or H**T is applied to C(1:m,1:n-k+i+ib-1) */

#line 331 "dormrq.f"
		ni = *n - *k + i__ + ib - 1;
#line 332 "dormrq.f"
	    }

/*           Apply H or H**T */

#line 336 "dormrq.f"
	    dlarfb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, &a[
		    i__ + a_dim1], lda, &work[iwt], &c__65, &c__[c_offset], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8, (
		    ftnlen)7);
#line 339 "dormrq.f"
/* L10: */
#line 339 "dormrq.f"
	}
#line 340 "dormrq.f"
    }
#line 341 "dormrq.f"
    work[1] = (doublereal) lwkopt;
#line 342 "dormrq.f"
    return 0;

/*     End of DORMRQ */

} /* dormrq_ */

