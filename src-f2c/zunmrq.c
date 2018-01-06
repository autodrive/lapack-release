#line 1 "zunmrq.f"
/* zunmrq.f -- translated by f2c (version 20100827).
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

#line 1 "zunmrq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b ZUNMRQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNMRQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmrq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmrq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmrq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > ZUNMRQ overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(1)**H H(2)**H . . . H(k)**H */
/* > */
/* > as returned by ZGERQF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          = 'C':  Transpose, apply Q**H. */
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
/* >          ZGERQF in the last k rows of its array argument A. */
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
/* >          reflector H(i), as returned by ZGERQF. */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zunmrq_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i1, i2, i3, ib, nb, mi, ni, nq, nw, iws;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int zunmr2_(char *, char *, integer *, integer *, 
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

#line 215 "zunmrq.f"
    /* Parameter adjustments */
#line 215 "zunmrq.f"
    a_dim1 = *lda;
#line 215 "zunmrq.f"
    a_offset = 1 + a_dim1;
#line 215 "zunmrq.f"
    a -= a_offset;
#line 215 "zunmrq.f"
    --tau;
#line 215 "zunmrq.f"
    c_dim1 = *ldc;
#line 215 "zunmrq.f"
    c_offset = 1 + c_dim1;
#line 215 "zunmrq.f"
    c__ -= c_offset;
#line 215 "zunmrq.f"
    --work;
#line 215 "zunmrq.f"

#line 215 "zunmrq.f"
    /* Function Body */
#line 215 "zunmrq.f"
    *info = 0;
#line 216 "zunmrq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 217 "zunmrq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 218 "zunmrq.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 222 "zunmrq.f"
    if (left) {
#line 223 "zunmrq.f"
	nq = *m;
#line 224 "zunmrq.f"
	nw = max(1,*n);
#line 225 "zunmrq.f"
    } else {
#line 226 "zunmrq.f"
	nq = *n;
#line 227 "zunmrq.f"
	nw = max(1,*m);
#line 228 "zunmrq.f"
    }
#line 229 "zunmrq.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 230 "zunmrq.f"
	*info = -1;
#line 231 "zunmrq.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 232 "zunmrq.f"
	*info = -2;
#line 233 "zunmrq.f"
    } else if (*m < 0) {
#line 234 "zunmrq.f"
	*info = -3;
#line 235 "zunmrq.f"
    } else if (*n < 0) {
#line 236 "zunmrq.f"
	*info = -4;
#line 237 "zunmrq.f"
    } else if (*k < 0 || *k > nq) {
#line 238 "zunmrq.f"
	*info = -5;
#line 239 "zunmrq.f"
    } else if (*lda < max(1,*k)) {
#line 240 "zunmrq.f"
	*info = -7;
#line 241 "zunmrq.f"
    } else if (*ldc < max(1,*m)) {
#line 242 "zunmrq.f"
	*info = -10;
#line 243 "zunmrq.f"
    }

#line 245 "zunmrq.f"
    if (*info == 0) {
#line 246 "zunmrq.f"
	if (*m == 0 || *n == 0) {
#line 247 "zunmrq.f"
	    lwkopt = 1;
#line 248 "zunmrq.f"
	} else {

/*           Determine the block size.  NB may be at most NBMAX, where */
/*           NBMAX is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 253 "zunmrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 253 "zunmrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 253 "zunmrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 253 "zunmrq.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 253 "zunmrq.f"
	    nb = min(i__1,i__2);
#line 255 "zunmrq.f"
	    lwkopt = nw * nb;
#line 256 "zunmrq.f"
	}
#line 257 "zunmrq.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 259 "zunmrq.f"
	if (*lwork < nw && ! lquery) {
#line 260 "zunmrq.f"
	    *info = -12;
#line 261 "zunmrq.f"
	}
#line 262 "zunmrq.f"
    }

#line 264 "zunmrq.f"
    if (*info != 0) {
#line 265 "zunmrq.f"
	i__1 = -(*info);
#line 265 "zunmrq.f"
	xerbla_("ZUNMRQ", &i__1, (ftnlen)6);
#line 266 "zunmrq.f"
	return 0;
#line 267 "zunmrq.f"
    } else if (lquery) {
#line 268 "zunmrq.f"
	return 0;
#line 269 "zunmrq.f"
    }

/*     Quick return if possible */

#line 273 "zunmrq.f"
    if (*m == 0 || *n == 0) {
#line 274 "zunmrq.f"
	return 0;
#line 275 "zunmrq.f"
    }

#line 277 "zunmrq.f"
    nbmin = 2;
#line 278 "zunmrq.f"
    ldwork = nw;
#line 279 "zunmrq.f"
    if (nb > 1 && nb < *k) {
#line 280 "zunmrq.f"
	iws = nw * nb;
#line 281 "zunmrq.f"
	if (*lwork < iws) {
#line 282 "zunmrq.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 283 "zunmrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 283 "zunmrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 283 "zunmrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 283 "zunmrq.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "ZUNMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 283 "zunmrq.f"
	    nbmin = max(i__1,i__2);
#line 285 "zunmrq.f"
	}
#line 286 "zunmrq.f"
    } else {
#line 287 "zunmrq.f"
	iws = nw;
#line 288 "zunmrq.f"
    }

#line 290 "zunmrq.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 294 "zunmrq.f"
	zunmr2_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 296 "zunmrq.f"
    } else {

/*        Use blocked code */

#line 300 "zunmrq.f"
	if (left && ! notran || ! left && notran) {
#line 302 "zunmrq.f"
	    i1 = 1;
#line 303 "zunmrq.f"
	    i2 = *k;
#line 304 "zunmrq.f"
	    i3 = nb;
#line 305 "zunmrq.f"
	} else {
#line 306 "zunmrq.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 307 "zunmrq.f"
	    i2 = 1;
#line 308 "zunmrq.f"
	    i3 = -nb;
#line 309 "zunmrq.f"
	}

#line 311 "zunmrq.f"
	if (left) {
#line 312 "zunmrq.f"
	    ni = *n;
#line 313 "zunmrq.f"
	} else {
#line 314 "zunmrq.f"
	    mi = *m;
#line 315 "zunmrq.f"
	}

#line 317 "zunmrq.f"
	if (notran) {
#line 318 "zunmrq.f"
	    *(unsigned char *)transt = 'C';
#line 319 "zunmrq.f"
	} else {
#line 320 "zunmrq.f"
	    *(unsigned char *)transt = 'N';
#line 321 "zunmrq.f"
	}

#line 323 "zunmrq.f"
	i__1 = i2;
#line 323 "zunmrq.f"
	i__2 = i3;
#line 323 "zunmrq.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 324 "zunmrq.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 324 "zunmrq.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 329 "zunmrq.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 329 "zunmrq.f"
	    zlarft_("Backward", "Rowwise", &i__4, &ib, &a[i__ + a_dim1], lda, 
		    &tau[i__], t, &c__65, (ftnlen)8, (ftnlen)7);
#line 331 "zunmrq.f"
	    if (left) {

/*              H or H**H is applied to C(1:m-k+i+ib-1,1:n) */

#line 335 "zunmrq.f"
		mi = *m - *k + i__ + ib - 1;
#line 336 "zunmrq.f"
	    } else {

/*              H or H**H is applied to C(1:m,1:n-k+i+ib-1) */

#line 340 "zunmrq.f"
		ni = *n - *k + i__ + ib - 1;
#line 341 "zunmrq.f"
	    }

/*           Apply H or H**H */

#line 345 "zunmrq.f"
	    zlarfb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, &a[
		    i__ + a_dim1], lda, t, &c__65, &c__[c_offset], ldc, &work[
		    1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8, (ftnlen)7);
#line 348 "zunmrq.f"
/* L10: */
#line 348 "zunmrq.f"
	}
#line 349 "zunmrq.f"
    }
#line 350 "zunmrq.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 351 "zunmrq.f"
    return 0;

/*     End of ZUNMRQ */

} /* zunmrq_ */

