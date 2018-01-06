#line 1 "sormql.f"
/* sormql.f -- translated by f2c (version 20100827).
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

#line 1 "sormql.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b SORMQL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORMQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormql.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormql.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormql.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > SORMQL overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by SGEQLF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is REAL array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          SGEQLF in the last k columns of its array argument A. */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sormql_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i1, i2, i3, ib, nb, mi, ni, nq, nw, iws;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int sorm2l_(char *, char *, integer *, integer *, 
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

#line 216 "sormql.f"
    /* Parameter adjustments */
#line 216 "sormql.f"
    a_dim1 = *lda;
#line 216 "sormql.f"
    a_offset = 1 + a_dim1;
#line 216 "sormql.f"
    a -= a_offset;
#line 216 "sormql.f"
    --tau;
#line 216 "sormql.f"
    c_dim1 = *ldc;
#line 216 "sormql.f"
    c_offset = 1 + c_dim1;
#line 216 "sormql.f"
    c__ -= c_offset;
#line 216 "sormql.f"
    --work;
#line 216 "sormql.f"

#line 216 "sormql.f"
    /* Function Body */
#line 216 "sormql.f"
    *info = 0;
#line 217 "sormql.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 218 "sormql.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 219 "sormql.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 223 "sormql.f"
    if (left) {
#line 224 "sormql.f"
	nq = *m;
#line 225 "sormql.f"
	nw = max(1,*n);
#line 226 "sormql.f"
    } else {
#line 227 "sormql.f"
	nq = *n;
#line 228 "sormql.f"
	nw = max(1,*m);
#line 229 "sormql.f"
    }
#line 230 "sormql.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 231 "sormql.f"
	*info = -1;
#line 232 "sormql.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 233 "sormql.f"
	*info = -2;
#line 234 "sormql.f"
    } else if (*m < 0) {
#line 235 "sormql.f"
	*info = -3;
#line 236 "sormql.f"
    } else if (*n < 0) {
#line 237 "sormql.f"
	*info = -4;
#line 238 "sormql.f"
    } else if (*k < 0 || *k > nq) {
#line 239 "sormql.f"
	*info = -5;
#line 240 "sormql.f"
    } else if (*lda < max(1,nq)) {
#line 241 "sormql.f"
	*info = -7;
#line 242 "sormql.f"
    } else if (*ldc < max(1,*m)) {
#line 243 "sormql.f"
	*info = -10;
#line 244 "sormql.f"
    }

#line 246 "sormql.f"
    if (*info == 0) {
#line 247 "sormql.f"
	if (*m == 0 || *n == 0) {
#line 248 "sormql.f"
	    lwkopt = 1;
#line 249 "sormql.f"
	} else {

/*           Determine the block size.  NB may be at most NBMAX, where */
/*           NBMAX is used to define the local array T. */


/* Computing MIN */
/* Writing concatenation */
#line 255 "sormql.f"
	    i__3[0] = 1, a__1[0] = side;
#line 255 "sormql.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 255 "sormql.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 255 "sormql.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "SORMQL", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 255 "sormql.f"
	    nb = min(i__1,i__2);
#line 257 "sormql.f"
	    lwkopt = nw * nb;
#line 258 "sormql.f"
	}
#line 259 "sormql.f"
	work[1] = (doublereal) lwkopt;

#line 261 "sormql.f"
	if (*lwork < nw && ! lquery) {
#line 262 "sormql.f"
	    *info = -12;
#line 263 "sormql.f"
	}
#line 264 "sormql.f"
    }

#line 266 "sormql.f"
    if (*info != 0) {
#line 267 "sormql.f"
	i__1 = -(*info);
#line 267 "sormql.f"
	xerbla_("SORMQL", &i__1, (ftnlen)6);
#line 268 "sormql.f"
	return 0;
#line 269 "sormql.f"
    } else if (lquery) {
#line 270 "sormql.f"
	return 0;
#line 271 "sormql.f"
    }

/*     Quick return if possible */

#line 275 "sormql.f"
    if (*m == 0 || *n == 0) {
#line 276 "sormql.f"
	return 0;
#line 277 "sormql.f"
    }

#line 279 "sormql.f"
    nbmin = 2;
#line 280 "sormql.f"
    ldwork = nw;
#line 281 "sormql.f"
    if (nb > 1 && nb < *k) {
#line 282 "sormql.f"
	iws = nw * nb;
#line 283 "sormql.f"
	if (*lwork < iws) {
#line 284 "sormql.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 285 "sormql.f"
	    i__3[0] = 1, a__1[0] = side;
#line 285 "sormql.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 285 "sormql.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 285 "sormql.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "SORMQL", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 285 "sormql.f"
	    nbmin = max(i__1,i__2);
#line 287 "sormql.f"
	}
#line 288 "sormql.f"
    } else {
#line 289 "sormql.f"
	iws = nw;
#line 290 "sormql.f"
    }

#line 292 "sormql.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 296 "sormql.f"
	sorm2l_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 298 "sormql.f"
    } else {

/*        Use blocked code */

#line 302 "sormql.f"
	if (left && notran || ! left && ! notran) {
#line 304 "sormql.f"
	    i1 = 1;
#line 305 "sormql.f"
	    i2 = *k;
#line 306 "sormql.f"
	    i3 = nb;
#line 307 "sormql.f"
	} else {
#line 308 "sormql.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 309 "sormql.f"
	    i2 = 1;
#line 310 "sormql.f"
	    i3 = -nb;
#line 311 "sormql.f"
	}

#line 313 "sormql.f"
	if (left) {
#line 314 "sormql.f"
	    ni = *n;
#line 315 "sormql.f"
	} else {
#line 316 "sormql.f"
	    mi = *m;
#line 317 "sormql.f"
	}

#line 319 "sormql.f"
	i__1 = i2;
#line 319 "sormql.f"
	i__2 = i3;
#line 319 "sormql.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 320 "sormql.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 320 "sormql.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 325 "sormql.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 325 "sormql.f"
	    slarft_("Backward", "Columnwise", &i__4, &ib, &a[i__ * a_dim1 + 1]
		    , lda, &tau[i__], t, &c__65, (ftnlen)8, (ftnlen)10);
#line 327 "sormql.f"
	    if (left) {

/*              H or H**T is applied to C(1:m-k+i+ib-1,1:n) */

#line 331 "sormql.f"
		mi = *m - *k + i__ + ib - 1;
#line 332 "sormql.f"
	    } else {

/*              H or H**T is applied to C(1:m,1:n-k+i+ib-1) */

#line 336 "sormql.f"
		ni = *n - *k + i__ + ib - 1;
#line 337 "sormql.f"
	    }

/*           Apply H or H**T */

#line 341 "sormql.f"
	    slarfb_(side, trans, "Backward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ * a_dim1 + 1], lda, t, &c__65, &c__[c_offset], ldc, &
		    work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8, (
		    ftnlen)10);
#line 344 "sormql.f"
/* L10: */
#line 344 "sormql.f"
	}
#line 345 "sormql.f"
    }
#line 346 "sormql.f"
    work[1] = (doublereal) lwkopt;
#line 347 "sormql.f"
    return 0;

/*     End of SORMQL */

} /* sormql_ */

