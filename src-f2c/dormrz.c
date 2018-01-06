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
    static integer i__;
    static doublereal t[4160]	/* was [65][64] */;
    static integer i1, i2, i3, ib, ic, ja, jc, nb, mi, ni, nq, nw, iws;
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

#line 235 "dormrz.f"
    /* Parameter adjustments */
#line 235 "dormrz.f"
    a_dim1 = *lda;
#line 235 "dormrz.f"
    a_offset = 1 + a_dim1;
#line 235 "dormrz.f"
    a -= a_offset;
#line 235 "dormrz.f"
    --tau;
#line 235 "dormrz.f"
    c_dim1 = *ldc;
#line 235 "dormrz.f"
    c_offset = 1 + c_dim1;
#line 235 "dormrz.f"
    c__ -= c_offset;
#line 235 "dormrz.f"
    --work;
#line 235 "dormrz.f"

#line 235 "dormrz.f"
    /* Function Body */
#line 235 "dormrz.f"
    *info = 0;
#line 236 "dormrz.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 237 "dormrz.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 238 "dormrz.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 242 "dormrz.f"
    if (left) {
#line 243 "dormrz.f"
	nq = *m;
#line 244 "dormrz.f"
	nw = max(1,*n);
#line 245 "dormrz.f"
    } else {
#line 246 "dormrz.f"
	nq = *n;
#line 247 "dormrz.f"
	nw = max(1,*m);
#line 248 "dormrz.f"
    }
#line 249 "dormrz.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 250 "dormrz.f"
	*info = -1;
#line 251 "dormrz.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 252 "dormrz.f"
	*info = -2;
#line 253 "dormrz.f"
    } else if (*m < 0) {
#line 254 "dormrz.f"
	*info = -3;
#line 255 "dormrz.f"
    } else if (*n < 0) {
#line 256 "dormrz.f"
	*info = -4;
#line 257 "dormrz.f"
    } else if (*k < 0 || *k > nq) {
#line 258 "dormrz.f"
	*info = -5;
#line 259 "dormrz.f"
    } else if (*l < 0 || left && *l > *m || ! left && *l > *n) {
#line 261 "dormrz.f"
	*info = -6;
#line 262 "dormrz.f"
    } else if (*lda < max(1,*k)) {
#line 263 "dormrz.f"
	*info = -8;
#line 264 "dormrz.f"
    } else if (*ldc < max(1,*m)) {
#line 265 "dormrz.f"
	*info = -11;
#line 266 "dormrz.f"
    }

#line 268 "dormrz.f"
    if (*info == 0) {
#line 269 "dormrz.f"
	if (*m == 0 || *n == 0) {
#line 270 "dormrz.f"
	    lwkopt = 1;
#line 271 "dormrz.f"
	} else {

/*           Determine the block size.  NB may be at most NBMAX, where */
/*           NBMAX is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 276 "dormrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 276 "dormrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 276 "dormrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 276 "dormrz.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 276 "dormrz.f"
	    nb = min(i__1,i__2);
#line 278 "dormrz.f"
	    lwkopt = nw * nb;
#line 279 "dormrz.f"
	}
#line 280 "dormrz.f"
	work[1] = (doublereal) lwkopt;

#line 282 "dormrz.f"
	if (*lwork < max(1,nw) && ! lquery) {
#line 283 "dormrz.f"
	    *info = -13;
#line 284 "dormrz.f"
	}
#line 285 "dormrz.f"
    }

#line 287 "dormrz.f"
    if (*info != 0) {
#line 288 "dormrz.f"
	i__1 = -(*info);
#line 288 "dormrz.f"
	xerbla_("DORMRZ", &i__1, (ftnlen)6);
#line 289 "dormrz.f"
	return 0;
#line 290 "dormrz.f"
    } else if (lquery) {
#line 291 "dormrz.f"
	return 0;
#line 292 "dormrz.f"
    }

/*     Quick return if possible */

#line 296 "dormrz.f"
    if (*m == 0 || *n == 0) {
#line 297 "dormrz.f"
	work[1] = 1.;
#line 298 "dormrz.f"
	return 0;
#line 299 "dormrz.f"
    }

#line 301 "dormrz.f"
    nbmin = 2;
#line 302 "dormrz.f"
    ldwork = nw;
#line 303 "dormrz.f"
    if (nb > 1 && nb < *k) {
#line 304 "dormrz.f"
	iws = nw * nb;
#line 305 "dormrz.f"
	if (*lwork < iws) {
#line 306 "dormrz.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 307 "dormrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 307 "dormrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 307 "dormrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 307 "dormrz.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 307 "dormrz.f"
	    nbmin = max(i__1,i__2);
#line 309 "dormrz.f"
	}
#line 310 "dormrz.f"
    } else {
#line 311 "dormrz.f"
	iws = nw;
#line 312 "dormrz.f"
    }

#line 314 "dormrz.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 318 "dormrz.f"
	dormr3_(side, trans, m, n, k, l, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 320 "dormrz.f"
    } else {

/*        Use blocked code */

#line 324 "dormrz.f"
	if (left && ! notran || ! left && notran) {
#line 326 "dormrz.f"
	    i1 = 1;
#line 327 "dormrz.f"
	    i2 = *k;
#line 328 "dormrz.f"
	    i3 = nb;
#line 329 "dormrz.f"
	} else {
#line 330 "dormrz.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 331 "dormrz.f"
	    i2 = 1;
#line 332 "dormrz.f"
	    i3 = -nb;
#line 333 "dormrz.f"
	}

#line 335 "dormrz.f"
	if (left) {
#line 336 "dormrz.f"
	    ni = *n;
#line 337 "dormrz.f"
	    jc = 1;
#line 338 "dormrz.f"
	    ja = *m - *l + 1;
#line 339 "dormrz.f"
	} else {
#line 340 "dormrz.f"
	    mi = *m;
#line 341 "dormrz.f"
	    ic = 1;
#line 342 "dormrz.f"
	    ja = *n - *l + 1;
#line 343 "dormrz.f"
	}

#line 345 "dormrz.f"
	if (notran) {
#line 346 "dormrz.f"
	    *(unsigned char *)transt = 'T';
#line 347 "dormrz.f"
	} else {
#line 348 "dormrz.f"
	    *(unsigned char *)transt = 'N';
#line 349 "dormrz.f"
	}

#line 351 "dormrz.f"
	i__1 = i2;
#line 351 "dormrz.f"
	i__2 = i3;
#line 351 "dormrz.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 352 "dormrz.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 352 "dormrz.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 357 "dormrz.f"
	    dlarzt_("Backward", "Rowwise", l, &ib, &a[i__ + ja * a_dim1], lda,
		     &tau[i__], t, &c__65, (ftnlen)8, (ftnlen)7);

#line 360 "dormrz.f"
	    if (left) {

/*              H or H**T is applied to C(i:m,1:n) */

#line 364 "dormrz.f"
		mi = *m - i__ + 1;
#line 365 "dormrz.f"
		ic = i__;
#line 366 "dormrz.f"
	    } else {

/*              H or H**T is applied to C(1:m,i:n) */

#line 370 "dormrz.f"
		ni = *n - i__ + 1;
#line 371 "dormrz.f"
		jc = i__;
#line 372 "dormrz.f"
	    }

/*           Apply H or H**T */

#line 376 "dormrz.f"
	    dlarzb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, l, &a[
		    i__ + ja * a_dim1], lda, t, &c__65, &c__[ic + jc * c_dim1]
		    , ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8,
		     (ftnlen)7);
#line 379 "dormrz.f"
/* L10: */
#line 379 "dormrz.f"
	}

#line 381 "dormrz.f"
    }

#line 383 "dormrz.f"
    work[1] = (doublereal) lwkopt;

#line 385 "dormrz.f"
    return 0;

/*     End of DORMRZ */

} /* dormrz_ */

