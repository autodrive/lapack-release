#line 1 "dormrz.f"
/* dormrz.f -- translated by f2c (version 20100827).
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

#line 1 "dormrz.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b DORMRZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORMRZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormrz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormrz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormrz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, L, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORMRZ overwrites the general real M-by-N matrix C with */
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
/* > as returned by DTZRZF. Q is of order M if SIDE = 'L' and of order N */
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
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The number of columns of the matrix A containing */
/* >          the meaningful part of the Householder reflectors. */
/* >          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DTZRZF in the last k rows of its array argument A. */
/* >          A is modified by the routine but restored on exit. */
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
/* >          reflector H(i), as returned by DTZRZF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
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

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dormrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
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
    static integer i__, i1, i2, i3, ib, ic, ja, jc, nb, mi, ni, nq, nw, iwt;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int dormr3_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen),
	     xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlarzb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlarzt_(char *, char 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
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

#line 231 "dormrz.f"
    /* Parameter adjustments */
#line 231 "dormrz.f"
    a_dim1 = *lda;
#line 231 "dormrz.f"
    a_offset = 1 + a_dim1;
#line 231 "dormrz.f"
    a -= a_offset;
#line 231 "dormrz.f"
    --tau;
#line 231 "dormrz.f"
    c_dim1 = *ldc;
#line 231 "dormrz.f"
    c_offset = 1 + c_dim1;
#line 231 "dormrz.f"
    c__ -= c_offset;
#line 231 "dormrz.f"
    --work;
#line 231 "dormrz.f"

#line 231 "dormrz.f"
    /* Function Body */
#line 231 "dormrz.f"
    *info = 0;
#line 232 "dormrz.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 233 "dormrz.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 234 "dormrz.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 238 "dormrz.f"
    if (left) {
#line 239 "dormrz.f"
	nq = *m;
#line 240 "dormrz.f"
	nw = max(1,*n);
#line 241 "dormrz.f"
    } else {
#line 242 "dormrz.f"
	nq = *n;
#line 243 "dormrz.f"
	nw = max(1,*m);
#line 244 "dormrz.f"
    }
#line 245 "dormrz.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 246 "dormrz.f"
	*info = -1;
#line 247 "dormrz.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 248 "dormrz.f"
	*info = -2;
#line 249 "dormrz.f"
    } else if (*m < 0) {
#line 250 "dormrz.f"
	*info = -3;
#line 251 "dormrz.f"
    } else if (*n < 0) {
#line 252 "dormrz.f"
	*info = -4;
#line 253 "dormrz.f"
    } else if (*k < 0 || *k > nq) {
#line 254 "dormrz.f"
	*info = -5;
#line 255 "dormrz.f"
    } else if (*l < 0 || left && *l > *m || ! left && *l > *n) {
#line 257 "dormrz.f"
	*info = -6;
#line 258 "dormrz.f"
    } else if (*lda < max(1,*k)) {
#line 259 "dormrz.f"
	*info = -8;
#line 260 "dormrz.f"
    } else if (*ldc < max(1,*m)) {
#line 261 "dormrz.f"
	*info = -11;
#line 262 "dormrz.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 263 "dormrz.f"
	*info = -13;
#line 264 "dormrz.f"
    }

#line 266 "dormrz.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 270 "dormrz.f"
	if (*m == 0 || *n == 0) {
#line 271 "dormrz.f"
	    lwkopt = 1;
#line 272 "dormrz.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 273 "dormrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 273 "dormrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 273 "dormrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 273 "dormrz.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 273 "dormrz.f"
	    nb = min(i__1,i__2);
#line 275 "dormrz.f"
	    lwkopt = nw * nb + 4160;
#line 276 "dormrz.f"
	}
#line 277 "dormrz.f"
	work[1] = (doublereal) lwkopt;
#line 278 "dormrz.f"
    }

#line 280 "dormrz.f"
    if (*info != 0) {
#line 281 "dormrz.f"
	i__1 = -(*info);
#line 281 "dormrz.f"
	xerbla_("DORMRZ", &i__1, (ftnlen)6);
#line 282 "dormrz.f"
	return 0;
#line 283 "dormrz.f"
    } else if (lquery) {
#line 284 "dormrz.f"
	return 0;
#line 285 "dormrz.f"
    }

/*     Quick return if possible */

#line 289 "dormrz.f"
    if (*m == 0 || *n == 0) {
#line 290 "dormrz.f"
	work[1] = 1.;
#line 291 "dormrz.f"
	return 0;
#line 292 "dormrz.f"
    }

#line 294 "dormrz.f"
    nbmin = 2;
#line 295 "dormrz.f"
    ldwork = nw;
#line 296 "dormrz.f"
    if (nb > 1 && nb < *k) {
#line 297 "dormrz.f"
	if (*lwork < nw * nb + 4160) {
#line 298 "dormrz.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 299 "dormrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 299 "dormrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 299 "dormrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 299 "dormrz.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 299 "dormrz.f"
	    nbmin = max(i__1,i__2);
#line 301 "dormrz.f"
	}
#line 302 "dormrz.f"
    }

#line 304 "dormrz.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 308 "dormrz.f"
	dormr3_(side, trans, m, n, k, l, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 310 "dormrz.f"
    } else {

/*        Use blocked code */

#line 314 "dormrz.f"
	iwt = nw * nb + 1;
#line 315 "dormrz.f"
	if (left && ! notran || ! left && notran) {
#line 317 "dormrz.f"
	    i1 = 1;
#line 318 "dormrz.f"
	    i2 = *k;
#line 319 "dormrz.f"
	    i3 = nb;
#line 320 "dormrz.f"
	} else {
#line 321 "dormrz.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 322 "dormrz.f"
	    i2 = 1;
#line 323 "dormrz.f"
	    i3 = -nb;
#line 324 "dormrz.f"
	}

#line 326 "dormrz.f"
	if (left) {
#line 327 "dormrz.f"
	    ni = *n;
#line 328 "dormrz.f"
	    jc = 1;
#line 329 "dormrz.f"
	    ja = *m - *l + 1;
#line 330 "dormrz.f"
	} else {
#line 331 "dormrz.f"
	    mi = *m;
#line 332 "dormrz.f"
	    ic = 1;
#line 333 "dormrz.f"
	    ja = *n - *l + 1;
#line 334 "dormrz.f"
	}

#line 336 "dormrz.f"
	if (notran) {
#line 337 "dormrz.f"
	    *(unsigned char *)transt = 'T';
#line 338 "dormrz.f"
	} else {
#line 339 "dormrz.f"
	    *(unsigned char *)transt = 'N';
#line 340 "dormrz.f"
	}

#line 342 "dormrz.f"
	i__1 = i2;
#line 342 "dormrz.f"
	i__2 = i3;
#line 342 "dormrz.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 343 "dormrz.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 343 "dormrz.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 348 "dormrz.f"
	    dlarzt_("Backward", "Rowwise", l, &ib, &a[i__ + ja * a_dim1], lda,
		     &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)7);

#line 351 "dormrz.f"
	    if (left) {

/*              H or H**T is applied to C(i:m,1:n) */

#line 355 "dormrz.f"
		mi = *m - i__ + 1;
#line 356 "dormrz.f"
		ic = i__;
#line 357 "dormrz.f"
	    } else {

/*              H or H**T is applied to C(1:m,i:n) */

#line 361 "dormrz.f"
		ni = *n - i__ + 1;
#line 362 "dormrz.f"
		jc = i__;
#line 363 "dormrz.f"
	    }

/*           Apply H or H**T */

#line 367 "dormrz.f"
	    dlarzb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, l, &a[
		    i__ + ja * a_dim1], lda, &work[iwt], &c__65, &c__[ic + jc 
		    * c_dim1], ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (
		    ftnlen)8, (ftnlen)7);
#line 370 "dormrz.f"
/* L10: */
#line 370 "dormrz.f"
	}

#line 372 "dormrz.f"
    }

#line 374 "dormrz.f"
    work[1] = (doublereal) lwkopt;

#line 376 "dormrz.f"
    return 0;

/*     End of DORMRZ */

} /* dormrz_ */

