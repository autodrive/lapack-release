#line 1 "sormrz.f"
/* sormrz.f -- translated by f2c (version 20100827).
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

#line 1 "sormrz.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b SORMRZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORMRZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormrz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormrz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormrz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, L, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORMRZ overwrites the general real M-by-N matrix C with */
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
/* > as returned by STZRZF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is REAL array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          STZRZF in the last k rows of its array argument A. */
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
/* >          TAU is REAL array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by STZRZF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
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
/* Subroutine */ int sormrz_(char *side, char *trans, integer *m, integer *n, 
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
    extern /* Subroutine */ int sormr3_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen),
	     xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarzb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical notran;
    static integer ldwork;
    static char transt[1];
    extern /* Subroutine */ int slarzt_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
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

#line 235 "sormrz.f"
    /* Parameter adjustments */
#line 235 "sormrz.f"
    a_dim1 = *lda;
#line 235 "sormrz.f"
    a_offset = 1 + a_dim1;
#line 235 "sormrz.f"
    a -= a_offset;
#line 235 "sormrz.f"
    --tau;
#line 235 "sormrz.f"
    c_dim1 = *ldc;
#line 235 "sormrz.f"
    c_offset = 1 + c_dim1;
#line 235 "sormrz.f"
    c__ -= c_offset;
#line 235 "sormrz.f"
    --work;
#line 235 "sormrz.f"

#line 235 "sormrz.f"
    /* Function Body */
#line 235 "sormrz.f"
    *info = 0;
#line 236 "sormrz.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 237 "sormrz.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 238 "sormrz.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 242 "sormrz.f"
    if (left) {
#line 243 "sormrz.f"
	nq = *m;
#line 244 "sormrz.f"
	nw = max(1,*n);
#line 245 "sormrz.f"
    } else {
#line 246 "sormrz.f"
	nq = *n;
#line 247 "sormrz.f"
	nw = max(1,*m);
#line 248 "sormrz.f"
    }
#line 249 "sormrz.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 250 "sormrz.f"
	*info = -1;
#line 251 "sormrz.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 252 "sormrz.f"
	*info = -2;
#line 253 "sormrz.f"
    } else if (*m < 0) {
#line 254 "sormrz.f"
	*info = -3;
#line 255 "sormrz.f"
    } else if (*n < 0) {
#line 256 "sormrz.f"
	*info = -4;
#line 257 "sormrz.f"
    } else if (*k < 0 || *k > nq) {
#line 258 "sormrz.f"
	*info = -5;
#line 259 "sormrz.f"
    } else if (*l < 0 || left && *l > *m || ! left && *l > *n) {
#line 261 "sormrz.f"
	*info = -6;
#line 262 "sormrz.f"
    } else if (*lda < max(1,*k)) {
#line 263 "sormrz.f"
	*info = -8;
#line 264 "sormrz.f"
    } else if (*ldc < max(1,*m)) {
#line 265 "sormrz.f"
	*info = -11;
#line 266 "sormrz.f"
    }

#line 268 "sormrz.f"
    if (*info == 0) {
#line 269 "sormrz.f"
	if (*m == 0 || *n == 0) {
#line 270 "sormrz.f"
	    lwkopt = 1;
#line 271 "sormrz.f"
	} else {

/*           Determine the block size.  NB may be at most NBMAX, where */
/*           NBMAX is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 276 "sormrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 276 "sormrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 276 "sormrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 276 "sormrz.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "SORMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 276 "sormrz.f"
	    nb = min(i__1,i__2);
#line 278 "sormrz.f"
	    lwkopt = nw * nb;
#line 279 "sormrz.f"
	}
#line 280 "sormrz.f"
	work[1] = (doublereal) lwkopt;

#line 282 "sormrz.f"
	if (*lwork < max(1,nw) && ! lquery) {
#line 283 "sormrz.f"
	    *info = -13;
#line 284 "sormrz.f"
	}
#line 285 "sormrz.f"
    }

#line 287 "sormrz.f"
    if (*info != 0) {
#line 288 "sormrz.f"
	i__1 = -(*info);
#line 288 "sormrz.f"
	xerbla_("SORMRZ", &i__1, (ftnlen)6);
#line 289 "sormrz.f"
	return 0;
#line 290 "sormrz.f"
    } else if (lquery) {
#line 291 "sormrz.f"
	return 0;
#line 292 "sormrz.f"
    }

/*     Quick return if possible */

#line 296 "sormrz.f"
    if (*m == 0 || *n == 0) {
#line 297 "sormrz.f"
	return 0;
#line 298 "sormrz.f"
    }

#line 300 "sormrz.f"
    nbmin = 2;
#line 301 "sormrz.f"
    ldwork = nw;
#line 302 "sormrz.f"
    if (nb > 1 && nb < *k) {
#line 303 "sormrz.f"
	iws = nw * nb;
#line 304 "sormrz.f"
	if (*lwork < iws) {
#line 305 "sormrz.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 306 "sormrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 306 "sormrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 306 "sormrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 306 "sormrz.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "SORMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 306 "sormrz.f"
	    nbmin = max(i__1,i__2);
#line 308 "sormrz.f"
	}
#line 309 "sormrz.f"
    } else {
#line 310 "sormrz.f"
	iws = nw;
#line 311 "sormrz.f"
    }

#line 313 "sormrz.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 317 "sormrz.f"
	sormr3_(side, trans, m, n, k, l, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 319 "sormrz.f"
    } else {

/*        Use blocked code */

#line 323 "sormrz.f"
	if (left && ! notran || ! left && notran) {
#line 325 "sormrz.f"
	    i1 = 1;
#line 326 "sormrz.f"
	    i2 = *k;
#line 327 "sormrz.f"
	    i3 = nb;
#line 328 "sormrz.f"
	} else {
#line 329 "sormrz.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 330 "sormrz.f"
	    i2 = 1;
#line 331 "sormrz.f"
	    i3 = -nb;
#line 332 "sormrz.f"
	}

#line 334 "sormrz.f"
	if (left) {
#line 335 "sormrz.f"
	    ni = *n;
#line 336 "sormrz.f"
	    jc = 1;
#line 337 "sormrz.f"
	    ja = *m - *l + 1;
#line 338 "sormrz.f"
	} else {
#line 339 "sormrz.f"
	    mi = *m;
#line 340 "sormrz.f"
	    ic = 1;
#line 341 "sormrz.f"
	    ja = *n - *l + 1;
#line 342 "sormrz.f"
	}

#line 344 "sormrz.f"
	if (notran) {
#line 345 "sormrz.f"
	    *(unsigned char *)transt = 'T';
#line 346 "sormrz.f"
	} else {
#line 347 "sormrz.f"
	    *(unsigned char *)transt = 'N';
#line 348 "sormrz.f"
	}

#line 350 "sormrz.f"
	i__1 = i2;
#line 350 "sormrz.f"
	i__2 = i3;
#line 350 "sormrz.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 351 "sormrz.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 351 "sormrz.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 356 "sormrz.f"
	    slarzt_("Backward", "Rowwise", l, &ib, &a[i__ + ja * a_dim1], lda,
		     &tau[i__], t, &c__65, (ftnlen)8, (ftnlen)7);

#line 359 "sormrz.f"
	    if (left) {

/*              H or H**T is applied to C(i:m,1:n) */

#line 363 "sormrz.f"
		mi = *m - i__ + 1;
#line 364 "sormrz.f"
		ic = i__;
#line 365 "sormrz.f"
	    } else {

/*              H or H**T is applied to C(1:m,i:n) */

#line 369 "sormrz.f"
		ni = *n - i__ + 1;
#line 370 "sormrz.f"
		jc = i__;
#line 371 "sormrz.f"
	    }

/*           Apply H or H**T */

#line 375 "sormrz.f"
	    slarzb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, l, &a[
		    i__ + ja * a_dim1], lda, t, &c__65, &c__[ic + jc * c_dim1]
		    , ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8,
		     (ftnlen)7);
#line 378 "sormrz.f"
/* L10: */
#line 378 "sormrz.f"
	}

#line 380 "sormrz.f"
    }

#line 382 "sormrz.f"
    work[1] = (doublereal) lwkopt;

#line 384 "sormrz.f"
    return 0;

/*     End of SORMRZ */

} /* sormrz_ */

