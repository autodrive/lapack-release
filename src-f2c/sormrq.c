#line 1 "sormrq.f"
/* sormrq.f -- translated by f2c (version 20100827).
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

#line 1 "sormrq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b SORMRQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORMRQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormrq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormrq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormrq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORMRQ overwrites the general real M-by-N matrix C with */
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
/* > as returned by SGERQF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is REAL array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          SGERQF in the last k rows of its array argument A. */
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
/* >          reflector H(i), as returned by SGERQF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sormrq_(char *side, char *trans, integer *m, integer *n, 
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
    extern /* Subroutine */ int sormr2_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen), slarfb_(char 
	    *, char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static logical notran;
    static integer ldwork;
    static char transt[1];
    static integer lwkopt;
    static logical lquery;


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

#line 213 "sormrq.f"
    /* Parameter adjustments */
#line 213 "sormrq.f"
    a_dim1 = *lda;
#line 213 "sormrq.f"
    a_offset = 1 + a_dim1;
#line 213 "sormrq.f"
    a -= a_offset;
#line 213 "sormrq.f"
    --tau;
#line 213 "sormrq.f"
    c_dim1 = *ldc;
#line 213 "sormrq.f"
    c_offset = 1 + c_dim1;
#line 213 "sormrq.f"
    c__ -= c_offset;
#line 213 "sormrq.f"
    --work;
#line 213 "sormrq.f"

#line 213 "sormrq.f"
    /* Function Body */
#line 213 "sormrq.f"
    *info = 0;
#line 214 "sormrq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 215 "sormrq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 216 "sormrq.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 220 "sormrq.f"
    if (left) {
#line 221 "sormrq.f"
	nq = *m;
#line 222 "sormrq.f"
	nw = max(1,*n);
#line 223 "sormrq.f"
    } else {
#line 224 "sormrq.f"
	nq = *n;
#line 225 "sormrq.f"
	nw = max(1,*m);
#line 226 "sormrq.f"
    }
#line 227 "sormrq.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 228 "sormrq.f"
	*info = -1;
#line 229 "sormrq.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 230 "sormrq.f"
	*info = -2;
#line 231 "sormrq.f"
    } else if (*m < 0) {
#line 232 "sormrq.f"
	*info = -3;
#line 233 "sormrq.f"
    } else if (*n < 0) {
#line 234 "sormrq.f"
	*info = -4;
#line 235 "sormrq.f"
    } else if (*k < 0 || *k > nq) {
#line 236 "sormrq.f"
	*info = -5;
#line 237 "sormrq.f"
    } else if (*lda < max(1,*k)) {
#line 238 "sormrq.f"
	*info = -7;
#line 239 "sormrq.f"
    } else if (*ldc < max(1,*m)) {
#line 240 "sormrq.f"
	*info = -10;
#line 241 "sormrq.f"
    } else if (*lwork < nw && ! lquery) {
#line 242 "sormrq.f"
	*info = -12;
#line 243 "sormrq.f"
    }

#line 245 "sormrq.f"
    if (*info == 0) {

/*     Compute the workspace requirements */

#line 249 "sormrq.f"
	if (*m == 0 || *n == 0) {
#line 250 "sormrq.f"
	    lwkopt = 1;
#line 251 "sormrq.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 252 "sormrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 252 "sormrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 252 "sormrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 252 "sormrq.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "SORMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 252 "sormrq.f"
	    nb = min(i__1,i__2);
#line 254 "sormrq.f"
	    lwkopt = nw * nb + 4160;
#line 255 "sormrq.f"
	}
#line 256 "sormrq.f"
	work[1] = (doublereal) lwkopt;
#line 257 "sormrq.f"
    }

#line 259 "sormrq.f"
    if (*info != 0) {
#line 260 "sormrq.f"
	i__1 = -(*info);
#line 260 "sormrq.f"
	xerbla_("SORMRQ", &i__1, (ftnlen)6);
#line 261 "sormrq.f"
	return 0;
#line 262 "sormrq.f"
    } else if (lquery) {
#line 263 "sormrq.f"
	return 0;
#line 264 "sormrq.f"
    }

/*     Quick return if possible */

#line 268 "sormrq.f"
    if (*m == 0 || *n == 0) {
#line 269 "sormrq.f"
	return 0;
#line 270 "sormrq.f"
    }

#line 272 "sormrq.f"
    nbmin = 2;
#line 273 "sormrq.f"
    ldwork = nw;
#line 274 "sormrq.f"
    if (nb > 1 && nb < *k) {
#line 275 "sormrq.f"
	if (*lwork < nw * nb + 4160) {
#line 276 "sormrq.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 277 "sormrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 277 "sormrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 277 "sormrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 277 "sormrq.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "SORMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 277 "sormrq.f"
	    nbmin = max(i__1,i__2);
#line 279 "sormrq.f"
	}
#line 280 "sormrq.f"
    }

#line 282 "sormrq.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 286 "sormrq.f"
	sormr2_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 288 "sormrq.f"
    } else {

/*        Use blocked code */

#line 292 "sormrq.f"
	iwt = nw * nb + 1;
#line 293 "sormrq.f"
	if (left && ! notran || ! left && notran) {
#line 295 "sormrq.f"
	    i1 = 1;
#line 296 "sormrq.f"
	    i2 = *k;
#line 297 "sormrq.f"
	    i3 = nb;
#line 298 "sormrq.f"
	} else {
#line 299 "sormrq.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 300 "sormrq.f"
	    i2 = 1;
#line 301 "sormrq.f"
	    i3 = -nb;
#line 302 "sormrq.f"
	}

#line 304 "sormrq.f"
	if (left) {
#line 305 "sormrq.f"
	    ni = *n;
#line 306 "sormrq.f"
	} else {
#line 307 "sormrq.f"
	    mi = *m;
#line 308 "sormrq.f"
	}

#line 310 "sormrq.f"
	if (notran) {
#line 311 "sormrq.f"
	    *(unsigned char *)transt = 'T';
#line 312 "sormrq.f"
	} else {
#line 313 "sormrq.f"
	    *(unsigned char *)transt = 'N';
#line 314 "sormrq.f"
	}

#line 316 "sormrq.f"
	i__1 = i2;
#line 316 "sormrq.f"
	i__2 = i3;
#line 316 "sormrq.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 317 "sormrq.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 317 "sormrq.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 322 "sormrq.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 322 "sormrq.f"
	    slarft_("Backward", "Rowwise", &i__4, &ib, &a[i__ + a_dim1], lda, 
		    &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)7);
#line 324 "sormrq.f"
	    if (left) {

/*              H or H**T is applied to C(1:m-k+i+ib-1,1:n) */

#line 328 "sormrq.f"
		mi = *m - *k + i__ + ib - 1;
#line 329 "sormrq.f"
	    } else {

/*              H or H**T is applied to C(1:m,1:n-k+i+ib-1) */

#line 333 "sormrq.f"
		ni = *n - *k + i__ + ib - 1;
#line 334 "sormrq.f"
	    }

/*           Apply H or H**T */

#line 338 "sormrq.f"
	    slarfb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, &a[
		    i__ + a_dim1], lda, &work[iwt], &c__65, &c__[c_offset], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8, (
		    ftnlen)7);
#line 341 "sormrq.f"
/* L10: */
#line 341 "sormrq.f"
	}
#line 342 "sormrq.f"
    }
#line 343 "sormrq.f"
    work[1] = (doublereal) lwkopt;
#line 344 "sormrq.f"
    return 0;

/*     End of SORMRQ */

} /* sormrq_ */

