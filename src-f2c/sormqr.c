#line 1 "sormqr.f"
/* sormqr.f -- translated by f2c (version 20100827).
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

#line 1 "sormqr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b SORMQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORMQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > SORMQR overwrites the general real M-by-N matrix C with */
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
/* > as returned by SGEQRF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          SGEQRF in the first k columns of its array argument A. */
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
/* >          reflector H(i), as returned by SGEQRF. */
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
/* Subroutine */ int sormqr_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i__, i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iwt;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int sorm2r_(char *, char *, integer *, integer *, 
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

#line 212 "sormqr.f"
    /* Parameter adjustments */
#line 212 "sormqr.f"
    a_dim1 = *lda;
#line 212 "sormqr.f"
    a_offset = 1 + a_dim1;
#line 212 "sormqr.f"
    a -= a_offset;
#line 212 "sormqr.f"
    --tau;
#line 212 "sormqr.f"
    c_dim1 = *ldc;
#line 212 "sormqr.f"
    c_offset = 1 + c_dim1;
#line 212 "sormqr.f"
    c__ -= c_offset;
#line 212 "sormqr.f"
    --work;
#line 212 "sormqr.f"

#line 212 "sormqr.f"
    /* Function Body */
#line 212 "sormqr.f"
    *info = 0;
#line 213 "sormqr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 214 "sormqr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 215 "sormqr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 219 "sormqr.f"
    if (left) {
#line 220 "sormqr.f"
	nq = *m;
#line 221 "sormqr.f"
	nw = *n;
#line 222 "sormqr.f"
    } else {
#line 223 "sormqr.f"
	nq = *n;
#line 224 "sormqr.f"
	nw = *m;
#line 225 "sormqr.f"
    }
#line 226 "sormqr.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 227 "sormqr.f"
	*info = -1;
#line 228 "sormqr.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 229 "sormqr.f"
	*info = -2;
#line 230 "sormqr.f"
    } else if (*m < 0) {
#line 231 "sormqr.f"
	*info = -3;
#line 232 "sormqr.f"
    } else if (*n < 0) {
#line 233 "sormqr.f"
	*info = -4;
#line 234 "sormqr.f"
    } else if (*k < 0 || *k > nq) {
#line 235 "sormqr.f"
	*info = -5;
#line 236 "sormqr.f"
    } else if (*lda < max(1,nq)) {
#line 237 "sormqr.f"
	*info = -7;
#line 238 "sormqr.f"
    } else if (*ldc < max(1,*m)) {
#line 239 "sormqr.f"
	*info = -10;
#line 240 "sormqr.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 241 "sormqr.f"
	*info = -12;
#line 242 "sormqr.f"
    }

#line 244 "sormqr.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

/* Computing MIN */
/* Writing concatenation */
#line 248 "sormqr.f"
	i__3[0] = 1, a__1[0] = side;
#line 248 "sormqr.f"
	i__3[1] = 1, a__1[1] = trans;
#line 248 "sormqr.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 248 "sormqr.f"
	i__1 = 64, i__2 = ilaenv_(&c__1, "SORMQR", ch__1, m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)2);
#line 248 "sormqr.f"
	nb = min(i__1,i__2);
#line 250 "sormqr.f"
	lwkopt = max(1,nw) * nb + 4160;
#line 251 "sormqr.f"
	work[1] = (doublereal) lwkopt;
#line 252 "sormqr.f"
    }

#line 254 "sormqr.f"
    if (*info != 0) {
#line 255 "sormqr.f"
	i__1 = -(*info);
#line 255 "sormqr.f"
	xerbla_("SORMQR", &i__1, (ftnlen)6);
#line 256 "sormqr.f"
	return 0;
#line 257 "sormqr.f"
    } else if (lquery) {
#line 258 "sormqr.f"
	return 0;
#line 259 "sormqr.f"
    }

/*     Quick return if possible */

#line 263 "sormqr.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 264 "sormqr.f"
	work[1] = 1.;
#line 265 "sormqr.f"
	return 0;
#line 266 "sormqr.f"
    }

#line 268 "sormqr.f"
    nbmin = 2;
#line 269 "sormqr.f"
    ldwork = nw;
#line 270 "sormqr.f"
    if (nb > 1 && nb < *k) {
#line 271 "sormqr.f"
	if (*lwork < nw * nb + 4160) {
#line 272 "sormqr.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 273 "sormqr.f"
	    i__3[0] = 1, a__1[0] = side;
#line 273 "sormqr.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 273 "sormqr.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 273 "sormqr.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "SORMQR", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 273 "sormqr.f"
	    nbmin = max(i__1,i__2);
#line 275 "sormqr.f"
	}
#line 276 "sormqr.f"
    }

#line 278 "sormqr.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 282 "sormqr.f"
	sorm2r_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 284 "sormqr.f"
    } else {

/*        Use blocked code */

#line 288 "sormqr.f"
	iwt = nw * nb + 1;
#line 289 "sormqr.f"
	if (left && ! notran || ! left && notran) {
#line 291 "sormqr.f"
	    i1 = 1;
#line 292 "sormqr.f"
	    i2 = *k;
#line 293 "sormqr.f"
	    i3 = nb;
#line 294 "sormqr.f"
	} else {
#line 295 "sormqr.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 296 "sormqr.f"
	    i2 = 1;
#line 297 "sormqr.f"
	    i3 = -nb;
#line 298 "sormqr.f"
	}

#line 300 "sormqr.f"
	if (left) {
#line 301 "sormqr.f"
	    ni = *n;
#line 302 "sormqr.f"
	    jc = 1;
#line 303 "sormqr.f"
	} else {
#line 304 "sormqr.f"
	    mi = *m;
#line 305 "sormqr.f"
	    ic = 1;
#line 306 "sormqr.f"
	}

#line 308 "sormqr.f"
	i__1 = i2;
#line 308 "sormqr.f"
	i__2 = i3;
#line 308 "sormqr.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 309 "sormqr.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 309 "sormqr.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i) H(i+1) . . . H(i+ib-1) */

#line 314 "sormqr.f"
	    i__4 = nq - i__ + 1;
#line 314 "sormqr.f"
	    slarft_("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &tau[i__], &work[iwt], &c__65, (ftnlen)7, (
		    ftnlen)10);
#line 316 "sormqr.f"
	    if (left) {

/*              H or H**T is applied to C(i:m,1:n) */

#line 320 "sormqr.f"
		mi = *m - i__ + 1;
#line 321 "sormqr.f"
		ic = i__;
#line 322 "sormqr.f"
	    } else {

/*              H or H**T is applied to C(1:m,i:n) */

#line 326 "sormqr.f"
		ni = *n - i__ + 1;
#line 327 "sormqr.f"
		jc = i__;
#line 328 "sormqr.f"
	    }

/*           Apply H or H**T */

#line 332 "sormqr.f"
	    slarfb_(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ + i__ * a_dim1], lda, &work[iwt], &c__65, &c__[ic + 
		    jc * c_dim1], ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)
		    1, (ftnlen)7, (ftnlen)10);
#line 335 "sormqr.f"
/* L10: */
#line 335 "sormqr.f"
	}
#line 336 "sormqr.f"
    }
#line 337 "sormqr.f"
    work[1] = (doublereal) lwkopt;
#line 338 "sormqr.f"
    return 0;

/*     End of SORMQR */

} /* sormqr_ */

