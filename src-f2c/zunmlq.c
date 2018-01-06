#line 1 "zunmlq.f"
/* zunmlq.f -- translated by f2c (version 20100827).
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

#line 1 "zunmlq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b ZUNMLQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNMLQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmlq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmlq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmlq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNMLQ overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(k)**H . . . H(2)**H H(1)**H */
/* > */
/* > as returned by ZGELQF. Q is of order M if SIDE = 'L' and of order N */
/* > if SIDE = 'R'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**H from the Left; */
/* >          = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'C':  Conjugate transpose, apply Q**H. */
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
/* >          A is COMPLEX*16 array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          ZGELQF in the first k rows of its array argument A. */
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
/* >          reflector H(i), as returned by ZGELQF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If SIDE = 'L', LWORK >= max(1,N); */
/* >          if SIDE = 'R', LWORK >= max(1,M). */
/* >          For optimum performance LWORK >= N*NB if SIDE 'L', and */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zunmlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	 integer *info, ftnlen side_len, ftnlen trans_len)
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
    static doublecomplex t[4160]	/* was [65][64] */;
    static integer i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int zunml2_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical notran;
    static integer ldwork;
    extern /* Subroutine */ int zlarft_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static char transt[1];
    static integer lwkopt;
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

#line 215 "zunmlq.f"
    /* Parameter adjustments */
#line 215 "zunmlq.f"
    a_dim1 = *lda;
#line 215 "zunmlq.f"
    a_offset = 1 + a_dim1;
#line 215 "zunmlq.f"
    a -= a_offset;
#line 215 "zunmlq.f"
    --tau;
#line 215 "zunmlq.f"
    c_dim1 = *ldc;
#line 215 "zunmlq.f"
    c_offset = 1 + c_dim1;
#line 215 "zunmlq.f"
    c__ -= c_offset;
#line 215 "zunmlq.f"
    --work;
#line 215 "zunmlq.f"

#line 215 "zunmlq.f"
    /* Function Body */
#line 215 "zunmlq.f"
    *info = 0;
#line 216 "zunmlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 217 "zunmlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 218 "zunmlq.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 222 "zunmlq.f"
    if (left) {
#line 223 "zunmlq.f"
	nq = *m;
#line 224 "zunmlq.f"
	nw = *n;
#line 225 "zunmlq.f"
    } else {
#line 226 "zunmlq.f"
	nq = *n;
#line 227 "zunmlq.f"
	nw = *m;
#line 228 "zunmlq.f"
    }
#line 229 "zunmlq.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 230 "zunmlq.f"
	*info = -1;
#line 231 "zunmlq.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 232 "zunmlq.f"
	*info = -2;
#line 233 "zunmlq.f"
    } else if (*m < 0) {
#line 234 "zunmlq.f"
	*info = -3;
#line 235 "zunmlq.f"
    } else if (*n < 0) {
#line 236 "zunmlq.f"
	*info = -4;
#line 237 "zunmlq.f"
    } else if (*k < 0 || *k > nq) {
#line 238 "zunmlq.f"
	*info = -5;
#line 239 "zunmlq.f"
    } else if (*lda < max(1,*k)) {
#line 240 "zunmlq.f"
	*info = -7;
#line 241 "zunmlq.f"
    } else if (*ldc < max(1,*m)) {
#line 242 "zunmlq.f"
	*info = -10;
#line 243 "zunmlq.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 244 "zunmlq.f"
	*info = -12;
#line 245 "zunmlq.f"
    }

#line 247 "zunmlq.f"
    if (*info == 0) {

/*        Determine the block size.  NB may be at most NBMAX, where NBMAX */
/*        is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 252 "zunmlq.f"
	i__3[0] = 1, a__1[0] = side;
#line 252 "zunmlq.f"
	i__3[1] = 1, a__1[1] = trans;
#line 252 "zunmlq.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 252 "zunmlq.f"
	i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMLQ", ch__1, m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)2);
#line 252 "zunmlq.f"
	nb = min(i__1,i__2);
#line 254 "zunmlq.f"
	lwkopt = max(1,nw) * nb;
#line 255 "zunmlq.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 256 "zunmlq.f"
    }

#line 258 "zunmlq.f"
    if (*info != 0) {
#line 259 "zunmlq.f"
	i__1 = -(*info);
#line 259 "zunmlq.f"
	xerbla_("ZUNMLQ", &i__1, (ftnlen)6);
#line 260 "zunmlq.f"
	return 0;
#line 261 "zunmlq.f"
    } else if (lquery) {
#line 262 "zunmlq.f"
	return 0;
#line 263 "zunmlq.f"
    }

/*     Quick return if possible */

#line 267 "zunmlq.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 268 "zunmlq.f"
	work[1].r = 1., work[1].i = 0.;
#line 269 "zunmlq.f"
	return 0;
#line 270 "zunmlq.f"
    }

#line 272 "zunmlq.f"
    nbmin = 2;
#line 273 "zunmlq.f"
    ldwork = nw;
#line 274 "zunmlq.f"
    if (nb > 1 && nb < *k) {
#line 275 "zunmlq.f"
	iws = nw * nb;
#line 276 "zunmlq.f"
	if (*lwork < iws) {
#line 277 "zunmlq.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 278 "zunmlq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 278 "zunmlq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 278 "zunmlq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 278 "zunmlq.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "ZUNMLQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 278 "zunmlq.f"
	    nbmin = max(i__1,i__2);
#line 280 "zunmlq.f"
	}
#line 281 "zunmlq.f"
    } else {
#line 282 "zunmlq.f"
	iws = nw;
#line 283 "zunmlq.f"
    }

#line 285 "zunmlq.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 289 "zunmlq.f"
	zunml2_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 291 "zunmlq.f"
    } else {

/*        Use blocked code */

#line 295 "zunmlq.f"
	if (left && notran || ! left && ! notran) {
#line 297 "zunmlq.f"
	    i1 = 1;
#line 298 "zunmlq.f"
	    i2 = *k;
#line 299 "zunmlq.f"
	    i3 = nb;
#line 300 "zunmlq.f"
	} else {
#line 301 "zunmlq.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 302 "zunmlq.f"
	    i2 = 1;
#line 303 "zunmlq.f"
	    i3 = -nb;
#line 304 "zunmlq.f"
	}

#line 306 "zunmlq.f"
	if (left) {
#line 307 "zunmlq.f"
	    ni = *n;
#line 308 "zunmlq.f"
	    jc = 1;
#line 309 "zunmlq.f"
	} else {
#line 310 "zunmlq.f"
	    mi = *m;
#line 311 "zunmlq.f"
	    ic = 1;
#line 312 "zunmlq.f"
	}

#line 314 "zunmlq.f"
	if (notran) {
#line 315 "zunmlq.f"
	    *(unsigned char *)transt = 'C';
#line 316 "zunmlq.f"
	} else {
#line 317 "zunmlq.f"
	    *(unsigned char *)transt = 'N';
#line 318 "zunmlq.f"
	}

#line 320 "zunmlq.f"
	i__1 = i2;
#line 320 "zunmlq.f"
	i__2 = i3;
#line 320 "zunmlq.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 321 "zunmlq.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 321 "zunmlq.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i) H(i+1) . . . H(i+ib-1) */

#line 326 "zunmlq.f"
	    i__4 = nq - i__ + 1;
#line 326 "zunmlq.f"
	    zlarft_("Forward", "Rowwise", &i__4, &ib, &a[i__ + i__ * a_dim1], 
		    lda, &tau[i__], t, &c__65, (ftnlen)7, (ftnlen)7);
#line 328 "zunmlq.f"
	    if (left) {

/*              H or H**H is applied to C(i:m,1:n) */

#line 332 "zunmlq.f"
		mi = *m - i__ + 1;
#line 333 "zunmlq.f"
		ic = i__;
#line 334 "zunmlq.f"
	    } else {

/*              H or H**H is applied to C(1:m,i:n) */

#line 338 "zunmlq.f"
		ni = *n - i__ + 1;
#line 339 "zunmlq.f"
		jc = i__;
#line 340 "zunmlq.f"
	    }

/*           Apply H or H**H */

#line 344 "zunmlq.f"
	    zlarfb_(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &a[i__ 
		    + i__ * a_dim1], lda, t, &c__65, &c__[ic + jc * c_dim1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)7, (
		    ftnlen)7);
#line 347 "zunmlq.f"
/* L10: */
#line 347 "zunmlq.f"
	}
#line 348 "zunmlq.f"
    }
#line 349 "zunmlq.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 350 "zunmlq.f"
    return 0;

/*     End of ZUNMLQ */

} /* zunmlq_ */

