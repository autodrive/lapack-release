#line 1 "dormqr.f"
/* dormqr.f -- translated by f2c (version 20100827).
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

#line 1 "dormqr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b DORMQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORMQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > DORMQR overwrites the general real M-by-N matrix C with */
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
/* > as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DGEQRF in the first k columns of its array argument A. */
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
/* >          reflector H(i), as returned by DGEQRF. */
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
/* >          For optimum performance LWORK >= N*NB if SIDE = 'L', and */
/* >          LWORK >= M*NB if SIDE = 'R', where NB is the optimal */
/* >          blocksize. */
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

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dormqr_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i__;
    static doublereal t[4160]	/* was [65][64] */;
    static integer i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int dorm2r_(char *, char *, integer *, integer *, 
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
    static integer ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 214 "dormqr.f"
    /* Parameter adjustments */
#line 214 "dormqr.f"
    a_dim1 = *lda;
#line 214 "dormqr.f"
    a_offset = 1 + a_dim1;
#line 214 "dormqr.f"
    a -= a_offset;
#line 214 "dormqr.f"
    --tau;
#line 214 "dormqr.f"
    c_dim1 = *ldc;
#line 214 "dormqr.f"
    c_offset = 1 + c_dim1;
#line 214 "dormqr.f"
    c__ -= c_offset;
#line 214 "dormqr.f"
    --work;
#line 214 "dormqr.f"

#line 214 "dormqr.f"
    /* Function Body */
#line 214 "dormqr.f"
    *info = 0;
#line 215 "dormqr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 216 "dormqr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 217 "dormqr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 221 "dormqr.f"
    if (left) {
#line 222 "dormqr.f"
	nq = *m;
#line 223 "dormqr.f"
	nw = *n;
#line 224 "dormqr.f"
    } else {
#line 225 "dormqr.f"
	nq = *n;
#line 226 "dormqr.f"
	nw = *m;
#line 227 "dormqr.f"
    }
#line 228 "dormqr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 229 "dormqr.f"
	*info = -1;
#line 230 "dormqr.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 231 "dormqr.f"
	*info = -2;
#line 232 "dormqr.f"
    } else if (*m < 0) {
#line 233 "dormqr.f"
	*info = -3;
#line 234 "dormqr.f"
    } else if (*n < 0) {
#line 235 "dormqr.f"
	*info = -4;
#line 236 "dormqr.f"
    } else if (*k < 0 || *k > nq) {
#line 237 "dormqr.f"
	*info = -5;
#line 238 "dormqr.f"
    } else if (*lda < max(1,nq)) {
#line 239 "dormqr.f"
	*info = -7;
#line 240 "dormqr.f"
    } else if (*ldc < max(1,*m)) {
#line 241 "dormqr.f"
	*info = -10;
#line 242 "dormqr.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 243 "dormqr.f"
	*info = -12;
#line 244 "dormqr.f"
    }

#line 246 "dormqr.f"
    if (*info == 0) {

/*        Determine the block size.  NB may be at most NBMAX, where NBMAX */
/*        is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 251 "dormqr.f"
	i__3[0] = 1, a__1[0] = side;
#line 251 "dormqr.f"
	i__3[1] = 1, a__1[1] = trans;
#line 251 "dormqr.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 251 "dormqr.f"
	i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", ch__1, m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)2);
#line 251 "dormqr.f"
	nb = min(i__1,i__2);
#line 253 "dormqr.f"
	lwkopt = max(1,nw) * nb;
#line 254 "dormqr.f"
	work[1] = (doublereal) lwkopt;
#line 255 "dormqr.f"
    }

#line 257 "dormqr.f"
    if (*info != 0) {
#line 258 "dormqr.f"
	i__1 = -(*info);
#line 258 "dormqr.f"
	xerbla_("DORMQR", &i__1, (ftnlen)6);
#line 259 "dormqr.f"
	return 0;
#line 260 "dormqr.f"
    } else if (lquery) {
#line 261 "dormqr.f"
	return 0;
#line 262 "dormqr.f"
    }

/*     Quick return if possible */

#line 266 "dormqr.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 267 "dormqr.f"
	work[1] = 1.;
#line 268 "dormqr.f"
	return 0;
#line 269 "dormqr.f"
    }

#line 271 "dormqr.f"
    nbmin = 2;
#line 272 "dormqr.f"
    ldwork = nw;
#line 273 "dormqr.f"
    if (nb > 1 && nb < *k) {
#line 274 "dormqr.f"
	iws = nw * nb;
#line 275 "dormqr.f"
	if (*lwork < iws) {
#line 276 "dormqr.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 277 "dormqr.f"
	    i__3[0] = 1, a__1[0] = side;
#line 277 "dormqr.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 277 "dormqr.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 277 "dormqr.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMQR", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 277 "dormqr.f"
	    nbmin = max(i__1,i__2);
#line 279 "dormqr.f"
	}
#line 280 "dormqr.f"
    } else {
#line 281 "dormqr.f"
	iws = nw;
#line 282 "dormqr.f"
    }

#line 284 "dormqr.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 288 "dormqr.f"
	dorm2r_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 290 "dormqr.f"
    } else {

/*        Use blocked code */

#line 294 "dormqr.f"
	if (left && ! notran || ! left && notran) {
#line 296 "dormqr.f"
	    i1 = 1;
#line 297 "dormqr.f"
	    i2 = *k;
#line 298 "dormqr.f"
	    i3 = nb;
#line 299 "dormqr.f"
	} else {
#line 300 "dormqr.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 301 "dormqr.f"
	    i2 = 1;
#line 302 "dormqr.f"
	    i3 = -nb;
#line 303 "dormqr.f"
	}

#line 305 "dormqr.f"
	if (left) {
#line 306 "dormqr.f"
	    ni = *n;
#line 307 "dormqr.f"
	    jc = 1;
#line 308 "dormqr.f"
	} else {
#line 309 "dormqr.f"
	    mi = *m;
#line 310 "dormqr.f"
	    ic = 1;
#line 311 "dormqr.f"
	}

#line 313 "dormqr.f"
	i__1 = i2;
#line 313 "dormqr.f"
	i__2 = i3;
#line 313 "dormqr.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 314 "dormqr.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 314 "dormqr.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i) H(i+1) . . . H(i+ib-1) */

#line 319 "dormqr.f"
	    i__4 = nq - i__ + 1;
#line 319 "dormqr.f"
	    dlarft_("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &tau[i__], t, &c__65, (ftnlen)7, (ftnlen)10)
		    ;
#line 321 "dormqr.f"
	    if (left) {

/*              H or H**T is applied to C(i:m,1:n) */

#line 325 "dormqr.f"
		mi = *m - i__ + 1;
#line 326 "dormqr.f"
		ic = i__;
#line 327 "dormqr.f"
	    } else {

/*              H or H**T is applied to C(1:m,i:n) */

#line 331 "dormqr.f"
		ni = *n - i__ + 1;
#line 332 "dormqr.f"
		jc = i__;
#line 333 "dormqr.f"
	    }

/*           Apply H or H**T */

#line 337 "dormqr.f"
	    dlarfb_(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ + i__ * a_dim1], lda, t, &c__65, &c__[ic + jc * 
		    c_dim1], ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (
		    ftnlen)7, (ftnlen)10);
#line 340 "dormqr.f"
/* L10: */
#line 340 "dormqr.f"
	}
#line 341 "dormqr.f"
    }
#line 342 "dormqr.f"
    work[1] = (doublereal) lwkopt;
#line 343 "dormqr.f"
    return 0;

/*     End of DORMQR */

} /* dormqr_ */

