#line 1 "sormlq.f"
/* sormlq.f -- translated by f2c (version 20100827).
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

#line 1 "sormlq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b SORMLQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORMLQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormlq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormlq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormlq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > SORMLQ overwrites the general real M-by-N matrix C with */
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
/* > as returned by SGELQF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          SGELQF in the first k rows of its array argument A. */
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
/* >          reflector H(i), as returned by SGELQF. */
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
/* Subroutine */ int sormlq_(char *side, char *trans, integer *m, integer *n, 
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
    extern /* Subroutine */ int sorml2_(char *, char *, integer *, integer *, 
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

#line 217 "sormlq.f"
    /* Parameter adjustments */
#line 217 "sormlq.f"
    a_dim1 = *lda;
#line 217 "sormlq.f"
    a_offset = 1 + a_dim1;
#line 217 "sormlq.f"
    a -= a_offset;
#line 217 "sormlq.f"
    --tau;
#line 217 "sormlq.f"
    c_dim1 = *ldc;
#line 217 "sormlq.f"
    c_offset = 1 + c_dim1;
#line 217 "sormlq.f"
    c__ -= c_offset;
#line 217 "sormlq.f"
    --work;
#line 217 "sormlq.f"

#line 217 "sormlq.f"
    /* Function Body */
#line 217 "sormlq.f"
    *info = 0;
#line 218 "sormlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 219 "sormlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 220 "sormlq.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 224 "sormlq.f"
    if (left) {
#line 225 "sormlq.f"
	nq = *m;
#line 226 "sormlq.f"
	nw = *n;
#line 227 "sormlq.f"
    } else {
#line 228 "sormlq.f"
	nq = *n;
#line 229 "sormlq.f"
	nw = *m;
#line 230 "sormlq.f"
    }
#line 231 "sormlq.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 232 "sormlq.f"
	*info = -1;
#line 233 "sormlq.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 234 "sormlq.f"
	*info = -2;
#line 235 "sormlq.f"
    } else if (*m < 0) {
#line 236 "sormlq.f"
	*info = -3;
#line 237 "sormlq.f"
    } else if (*n < 0) {
#line 238 "sormlq.f"
	*info = -4;
#line 239 "sormlq.f"
    } else if (*k < 0 || *k > nq) {
#line 240 "sormlq.f"
	*info = -5;
#line 241 "sormlq.f"
    } else if (*lda < max(1,*k)) {
#line 242 "sormlq.f"
	*info = -7;
#line 243 "sormlq.f"
    } else if (*ldc < max(1,*m)) {
#line 244 "sormlq.f"
	*info = -10;
#line 245 "sormlq.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 246 "sormlq.f"
	*info = -12;
#line 247 "sormlq.f"
    }

#line 249 "sormlq.f"
    if (*info == 0) {

/*        Determine the block size.  NB may be at most NBMAX, where NBMAX */
/*        is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 254 "sormlq.f"
	i__3[0] = 1, a__1[0] = side;
#line 254 "sormlq.f"
	i__3[1] = 1, a__1[1] = trans;
#line 254 "sormlq.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 254 "sormlq.f"
	i__1 = 64, i__2 = ilaenv_(&c__1, "SORMLQ", ch__1, m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)2);
#line 254 "sormlq.f"
	nb = min(i__1,i__2);
#line 256 "sormlq.f"
	lwkopt = max(1,nw) * nb;
#line 257 "sormlq.f"
	work[1] = (doublereal) lwkopt;
#line 258 "sormlq.f"
    }

#line 260 "sormlq.f"
    if (*info != 0) {
#line 261 "sormlq.f"
	i__1 = -(*info);
#line 261 "sormlq.f"
	xerbla_("SORMLQ", &i__1, (ftnlen)6);
#line 262 "sormlq.f"
	return 0;
#line 263 "sormlq.f"
    } else if (lquery) {
#line 264 "sormlq.f"
	return 0;
#line 265 "sormlq.f"
    }

/*     Quick return if possible */

#line 269 "sormlq.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 270 "sormlq.f"
	work[1] = 1.;
#line 271 "sormlq.f"
	return 0;
#line 272 "sormlq.f"
    }

#line 274 "sormlq.f"
    nbmin = 2;
#line 275 "sormlq.f"
    ldwork = nw;
#line 276 "sormlq.f"
    if (nb > 1 && nb < *k) {
#line 277 "sormlq.f"
	iws = nw * nb;
#line 278 "sormlq.f"
	if (*lwork < iws) {
#line 279 "sormlq.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 280 "sormlq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 280 "sormlq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 280 "sormlq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 280 "sormlq.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "SORMLQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 280 "sormlq.f"
	    nbmin = max(i__1,i__2);
#line 282 "sormlq.f"
	}
#line 283 "sormlq.f"
    } else {
#line 284 "sormlq.f"
	iws = nw;
#line 285 "sormlq.f"
    }

#line 287 "sormlq.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 291 "sormlq.f"
	sorml2_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 293 "sormlq.f"
    } else {

/*        Use blocked code */

#line 297 "sormlq.f"
	if (left && notran || ! left && ! notran) {
#line 299 "sormlq.f"
	    i1 = 1;
#line 300 "sormlq.f"
	    i2 = *k;
#line 301 "sormlq.f"
	    i3 = nb;
#line 302 "sormlq.f"
	} else {
#line 303 "sormlq.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 304 "sormlq.f"
	    i2 = 1;
#line 305 "sormlq.f"
	    i3 = -nb;
#line 306 "sormlq.f"
	}

#line 308 "sormlq.f"
	if (left) {
#line 309 "sormlq.f"
	    ni = *n;
#line 310 "sormlq.f"
	    jc = 1;
#line 311 "sormlq.f"
	} else {
#line 312 "sormlq.f"
	    mi = *m;
#line 313 "sormlq.f"
	    ic = 1;
#line 314 "sormlq.f"
	}

#line 316 "sormlq.f"
	if (notran) {
#line 317 "sormlq.f"
	    *(unsigned char *)transt = 'T';
#line 318 "sormlq.f"
	} else {
#line 319 "sormlq.f"
	    *(unsigned char *)transt = 'N';
#line 320 "sormlq.f"
	}

#line 322 "sormlq.f"
	i__1 = i2;
#line 322 "sormlq.f"
	i__2 = i3;
#line 322 "sormlq.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 323 "sormlq.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 323 "sormlq.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i) H(i+1) . . . H(i+ib-1) */

#line 328 "sormlq.f"
	    i__4 = nq - i__ + 1;
#line 328 "sormlq.f"
	    slarft_("Forward", "Rowwise", &i__4, &ib, &a[i__ + i__ * a_dim1], 
		    lda, &tau[i__], t, &c__65, (ftnlen)7, (ftnlen)7);
#line 330 "sormlq.f"
	    if (left) {

/*              H or H**T is applied to C(i:m,1:n) */

#line 334 "sormlq.f"
		mi = *m - i__ + 1;
#line 335 "sormlq.f"
		ic = i__;
#line 336 "sormlq.f"
	    } else {

/*              H or H**T is applied to C(1:m,i:n) */

#line 340 "sormlq.f"
		ni = *n - i__ + 1;
#line 341 "sormlq.f"
		jc = i__;
#line 342 "sormlq.f"
	    }

/*           Apply H or H**T */

#line 346 "sormlq.f"
	    slarfb_(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &a[i__ 
		    + i__ * a_dim1], lda, t, &c__65, &c__[ic + jc * c_dim1], 
		    ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)7, (
		    ftnlen)7);
#line 349 "sormlq.f"
/* L10: */
#line 349 "sormlq.f"
	}
#line 350 "sormlq.f"
    }
#line 351 "sormlq.f"
    work[1] = (doublereal) lwkopt;
#line 352 "sormlq.f"
    return 0;

/*     End of SORMLQ */

} /* sormlq_ */

