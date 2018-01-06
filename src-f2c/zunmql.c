#line 1 "zunmql.f"
/* zunmql.f -- translated by f2c (version 20100827).
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

#line 1 "zunmql.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b ZUNMQL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNMQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmql.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmql.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmql.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > ZUNMQL overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by ZGEQLF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is COMPLEX*16 array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          ZGEQLF in the last k columns of its array argument A. */
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
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZGEQLF. */
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
/* Subroutine */ int zunmql_(char *side, char *trans, integer *m, integer *n, 
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
    extern /* Subroutine */ int zunm2l_(char *, char *, integer *, integer *, 
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

#line 214 "zunmql.f"
    /* Parameter adjustments */
#line 214 "zunmql.f"
    a_dim1 = *lda;
#line 214 "zunmql.f"
    a_offset = 1 + a_dim1;
#line 214 "zunmql.f"
    a -= a_offset;
#line 214 "zunmql.f"
    --tau;
#line 214 "zunmql.f"
    c_dim1 = *ldc;
#line 214 "zunmql.f"
    c_offset = 1 + c_dim1;
#line 214 "zunmql.f"
    c__ -= c_offset;
#line 214 "zunmql.f"
    --work;
#line 214 "zunmql.f"

#line 214 "zunmql.f"
    /* Function Body */
#line 214 "zunmql.f"
    *info = 0;
#line 215 "zunmql.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 216 "zunmql.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 217 "zunmql.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 221 "zunmql.f"
    if (left) {
#line 222 "zunmql.f"
	nq = *m;
#line 223 "zunmql.f"
	nw = max(1,*n);
#line 224 "zunmql.f"
    } else {
#line 225 "zunmql.f"
	nq = *n;
#line 226 "zunmql.f"
	nw = max(1,*m);
#line 227 "zunmql.f"
    }
#line 228 "zunmql.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 229 "zunmql.f"
	*info = -1;
#line 230 "zunmql.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 231 "zunmql.f"
	*info = -2;
#line 232 "zunmql.f"
    } else if (*m < 0) {
#line 233 "zunmql.f"
	*info = -3;
#line 234 "zunmql.f"
    } else if (*n < 0) {
#line 235 "zunmql.f"
	*info = -4;
#line 236 "zunmql.f"
    } else if (*k < 0 || *k > nq) {
#line 237 "zunmql.f"
	*info = -5;
#line 238 "zunmql.f"
    } else if (*lda < max(1,nq)) {
#line 239 "zunmql.f"
	*info = -7;
#line 240 "zunmql.f"
    } else if (*ldc < max(1,*m)) {
#line 241 "zunmql.f"
	*info = -10;
#line 242 "zunmql.f"
    }

#line 244 "zunmql.f"
    if (*info == 0) {
#line 245 "zunmql.f"
	if (*m == 0 || *n == 0) {
#line 246 "zunmql.f"
	    lwkopt = 1;
#line 247 "zunmql.f"
	} else {

/*           Determine the block size.  NB may be at most NBMAX, where */
/*           NBMAX is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 252 "zunmql.f"
	    i__3[0] = 1, a__1[0] = side;
#line 252 "zunmql.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 252 "zunmql.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 252 "zunmql.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMQL", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 252 "zunmql.f"
	    nb = min(i__1,i__2);
#line 254 "zunmql.f"
	    lwkopt = nw * nb;
#line 255 "zunmql.f"
	}
#line 256 "zunmql.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 258 "zunmql.f"
	if (*lwork < nw && ! lquery) {
#line 259 "zunmql.f"
	    *info = -12;
#line 260 "zunmql.f"
	}
#line 261 "zunmql.f"
    }

#line 263 "zunmql.f"
    if (*info != 0) {
#line 264 "zunmql.f"
	i__1 = -(*info);
#line 264 "zunmql.f"
	xerbla_("ZUNMQL", &i__1, (ftnlen)6);
#line 265 "zunmql.f"
	return 0;
#line 266 "zunmql.f"
    } else if (lquery) {
#line 267 "zunmql.f"
	return 0;
#line 268 "zunmql.f"
    }

/*     Quick return if possible */

#line 272 "zunmql.f"
    if (*m == 0 || *n == 0) {
#line 273 "zunmql.f"
	return 0;
#line 274 "zunmql.f"
    }

#line 276 "zunmql.f"
    nbmin = 2;
#line 277 "zunmql.f"
    ldwork = nw;
#line 278 "zunmql.f"
    if (nb > 1 && nb < *k) {
#line 279 "zunmql.f"
	iws = nw * nb;
#line 280 "zunmql.f"
	if (*lwork < iws) {
#line 281 "zunmql.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 282 "zunmql.f"
	    i__3[0] = 1, a__1[0] = side;
#line 282 "zunmql.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 282 "zunmql.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 282 "zunmql.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "ZUNMQL", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 282 "zunmql.f"
	    nbmin = max(i__1,i__2);
#line 284 "zunmql.f"
	}
#line 285 "zunmql.f"
    } else {
#line 286 "zunmql.f"
	iws = nw;
#line 287 "zunmql.f"
    }

#line 289 "zunmql.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 293 "zunmql.f"
	zunm2l_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 295 "zunmql.f"
    } else {

/*        Use blocked code */

#line 299 "zunmql.f"
	if (left && notran || ! left && ! notran) {
#line 301 "zunmql.f"
	    i1 = 1;
#line 302 "zunmql.f"
	    i2 = *k;
#line 303 "zunmql.f"
	    i3 = nb;
#line 304 "zunmql.f"
	} else {
#line 305 "zunmql.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 306 "zunmql.f"
	    i2 = 1;
#line 307 "zunmql.f"
	    i3 = -nb;
#line 308 "zunmql.f"
	}

#line 310 "zunmql.f"
	if (left) {
#line 311 "zunmql.f"
	    ni = *n;
#line 312 "zunmql.f"
	} else {
#line 313 "zunmql.f"
	    mi = *m;
#line 314 "zunmql.f"
	}

#line 316 "zunmql.f"
	i__1 = i2;
#line 316 "zunmql.f"
	i__2 = i3;
#line 316 "zunmql.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 317 "zunmql.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 317 "zunmql.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 322 "zunmql.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 322 "zunmql.f"
	    zlarft_("Backward", "Columnwise", &i__4, &ib, &a[i__ * a_dim1 + 1]
		    , lda, &tau[i__], t, &c__65, (ftnlen)8, (ftnlen)10);
#line 324 "zunmql.f"
	    if (left) {

/*              H or H**H is applied to C(1:m-k+i+ib-1,1:n) */

#line 328 "zunmql.f"
		mi = *m - *k + i__ + ib - 1;
#line 329 "zunmql.f"
	    } else {

/*              H or H**H is applied to C(1:m,1:n-k+i+ib-1) */

#line 333 "zunmql.f"
		ni = *n - *k + i__ + ib - 1;
#line 334 "zunmql.f"
	    }

/*           Apply H or H**H */

#line 338 "zunmql.f"
	    zlarfb_(side, trans, "Backward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ * a_dim1 + 1], lda, t, &c__65, &c__[c_offset], ldc, &
		    work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8, (
		    ftnlen)10);
#line 341 "zunmql.f"
/* L10: */
#line 341 "zunmql.f"
	}
#line 342 "zunmql.f"
    }
#line 343 "zunmql.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 344 "zunmql.f"
    return 0;

/*     End of ZUNMQL */

} /* zunmql_ */

