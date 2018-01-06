#line 1 "cunmrq.f"
/* cunmrq.f -- translated by f2c (version 20100827).
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

#line 1 "cunmrq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b CUNMRQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMRQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmrq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmrq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmrq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNMRQ overwrites the general complex M-by-N matrix C with */
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
/* > as returned by CGERQF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGERQF in the last k rows of its array argument A. */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunmrq_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i__, i1, i2, i3, ib, nb, mi, ni, nq, nw, iwt;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int cunmr2_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), clarfb_(char *, char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublecomplex *, integer *
	    , ftnlen, ftnlen, ftnlen, ftnlen), clarft_(char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
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

#line 213 "cunmrq.f"
    /* Parameter adjustments */
#line 213 "cunmrq.f"
    a_dim1 = *lda;
#line 213 "cunmrq.f"
    a_offset = 1 + a_dim1;
#line 213 "cunmrq.f"
    a -= a_offset;
#line 213 "cunmrq.f"
    --tau;
#line 213 "cunmrq.f"
    c_dim1 = *ldc;
#line 213 "cunmrq.f"
    c_offset = 1 + c_dim1;
#line 213 "cunmrq.f"
    c__ -= c_offset;
#line 213 "cunmrq.f"
    --work;
#line 213 "cunmrq.f"

#line 213 "cunmrq.f"
    /* Function Body */
#line 213 "cunmrq.f"
    *info = 0;
#line 214 "cunmrq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 215 "cunmrq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 216 "cunmrq.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 220 "cunmrq.f"
    if (left) {
#line 221 "cunmrq.f"
	nq = *m;
#line 222 "cunmrq.f"
	nw = max(1,*n);
#line 223 "cunmrq.f"
    } else {
#line 224 "cunmrq.f"
	nq = *n;
#line 225 "cunmrq.f"
	nw = max(1,*m);
#line 226 "cunmrq.f"
    }
#line 227 "cunmrq.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 228 "cunmrq.f"
	*info = -1;
#line 229 "cunmrq.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 230 "cunmrq.f"
	*info = -2;
#line 231 "cunmrq.f"
    } else if (*m < 0) {
#line 232 "cunmrq.f"
	*info = -3;
#line 233 "cunmrq.f"
    } else if (*n < 0) {
#line 234 "cunmrq.f"
	*info = -4;
#line 235 "cunmrq.f"
    } else if (*k < 0 || *k > nq) {
#line 236 "cunmrq.f"
	*info = -5;
#line 237 "cunmrq.f"
    } else if (*lda < max(1,*k)) {
#line 238 "cunmrq.f"
	*info = -7;
#line 239 "cunmrq.f"
    } else if (*ldc < max(1,*m)) {
#line 240 "cunmrq.f"
	*info = -10;
#line 241 "cunmrq.f"
    } else if (*lwork < nw && ! lquery) {
#line 242 "cunmrq.f"
	*info = -12;
#line 243 "cunmrq.f"
    }

#line 245 "cunmrq.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 249 "cunmrq.f"
	if (*m == 0 || *n == 0) {
#line 250 "cunmrq.f"
	    lwkopt = 1;
#line 251 "cunmrq.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 252 "cunmrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 252 "cunmrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 252 "cunmrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 252 "cunmrq.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "CUNMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 252 "cunmrq.f"
	    nb = min(i__1,i__2);
#line 254 "cunmrq.f"
	    lwkopt = nw * nb + 4160;
#line 255 "cunmrq.f"
	}
#line 256 "cunmrq.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 257 "cunmrq.f"
    }

#line 259 "cunmrq.f"
    if (*info != 0) {
#line 260 "cunmrq.f"
	i__1 = -(*info);
#line 260 "cunmrq.f"
	xerbla_("CUNMRQ", &i__1, (ftnlen)6);
#line 261 "cunmrq.f"
	return 0;
#line 262 "cunmrq.f"
    } else if (lquery) {
#line 263 "cunmrq.f"
	return 0;
#line 264 "cunmrq.f"
    }

/*     Quick return if possible */

#line 268 "cunmrq.f"
    if (*m == 0 || *n == 0) {
#line 269 "cunmrq.f"
	return 0;
#line 270 "cunmrq.f"
    }

#line 272 "cunmrq.f"
    nbmin = 2;
#line 273 "cunmrq.f"
    ldwork = nw;
#line 274 "cunmrq.f"
    if (nb > 1 && nb < *k) {
#line 275 "cunmrq.f"
	if (*lwork < nw * nb + 4160) {
#line 276 "cunmrq.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 277 "cunmrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 277 "cunmrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 277 "cunmrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 277 "cunmrq.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CUNMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 277 "cunmrq.f"
	    nbmin = max(i__1,i__2);
#line 279 "cunmrq.f"
	}
#line 280 "cunmrq.f"
    }

#line 282 "cunmrq.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 286 "cunmrq.f"
	cunmr2_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 288 "cunmrq.f"
    } else {

/*        Use blocked code */

#line 292 "cunmrq.f"
	iwt = nw * nb + 1;
#line 293 "cunmrq.f"
	if (left && ! notran || ! left && notran) {
#line 295 "cunmrq.f"
	    i1 = 1;
#line 296 "cunmrq.f"
	    i2 = *k;
#line 297 "cunmrq.f"
	    i3 = nb;
#line 298 "cunmrq.f"
	} else {
#line 299 "cunmrq.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 300 "cunmrq.f"
	    i2 = 1;
#line 301 "cunmrq.f"
	    i3 = -nb;
#line 302 "cunmrq.f"
	}

#line 304 "cunmrq.f"
	if (left) {
#line 305 "cunmrq.f"
	    ni = *n;
#line 306 "cunmrq.f"
	} else {
#line 307 "cunmrq.f"
	    mi = *m;
#line 308 "cunmrq.f"
	}

#line 310 "cunmrq.f"
	if (notran) {
#line 311 "cunmrq.f"
	    *(unsigned char *)transt = 'C';
#line 312 "cunmrq.f"
	} else {
#line 313 "cunmrq.f"
	    *(unsigned char *)transt = 'N';
#line 314 "cunmrq.f"
	}

#line 316 "cunmrq.f"
	i__1 = i2;
#line 316 "cunmrq.f"
	i__2 = i3;
#line 316 "cunmrq.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 317 "cunmrq.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 317 "cunmrq.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 322 "cunmrq.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 322 "cunmrq.f"
	    clarft_("Backward", "Rowwise", &i__4, &ib, &a[i__ + a_dim1], lda, 
		    &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)7);
#line 324 "cunmrq.f"
	    if (left) {

/*              H or H**H is applied to C(1:m-k+i+ib-1,1:n) */

#line 328 "cunmrq.f"
		mi = *m - *k + i__ + ib - 1;
#line 329 "cunmrq.f"
	    } else {

/*              H or H**H is applied to C(1:m,1:n-k+i+ib-1) */

#line 333 "cunmrq.f"
		ni = *n - *k + i__ + ib - 1;
#line 334 "cunmrq.f"
	    }

/*           Apply H or H**H */

#line 338 "cunmrq.f"
	    clarfb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, &a[
		    i__ + a_dim1], lda, &work[iwt], &c__65, &c__[c_offset], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8, (
		    ftnlen)7);
#line 341 "cunmrq.f"
/* L10: */
#line 341 "cunmrq.f"
	}
#line 342 "cunmrq.f"
    }
#line 343 "cunmrq.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 344 "cunmrq.f"
    return 0;

/*     End of CUNMRQ */

} /* cunmrq_ */

