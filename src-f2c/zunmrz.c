#line 1 "zunmrz.f"
/* zunmrz.f -- translated by f2c (version 20100827).
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

#line 1 "zunmrz.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b ZUNMRZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNMRZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmrz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmrz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmrz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, L, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNMRZ overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by ZTZRZF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          = 'C':  Conjugate transpose, apply Q**H. */
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
/* >          A is COMPLEX*16 array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          ZTZRZF in the last k rows of its array argument A. */
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
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZTZRZF. */
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
/* Subroutine */ int zunmrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	lwork, integer *info, ftnlen side_len, ftnlen trans_len)
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
    static integer i1, i2, i3, ib, ic, ja, jc, nb, mi, ni, nq, nw, iws;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int zunmr3_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical notran;
    static integer ldwork;
    extern /* Subroutine */ int zlarzb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static char transt[1];
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zlarzt_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);


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

#line 235 "zunmrz.f"
    /* Parameter adjustments */
#line 235 "zunmrz.f"
    a_dim1 = *lda;
#line 235 "zunmrz.f"
    a_offset = 1 + a_dim1;
#line 235 "zunmrz.f"
    a -= a_offset;
#line 235 "zunmrz.f"
    --tau;
#line 235 "zunmrz.f"
    c_dim1 = *ldc;
#line 235 "zunmrz.f"
    c_offset = 1 + c_dim1;
#line 235 "zunmrz.f"
    c__ -= c_offset;
#line 235 "zunmrz.f"
    --work;
#line 235 "zunmrz.f"

#line 235 "zunmrz.f"
    /* Function Body */
#line 235 "zunmrz.f"
    *info = 0;
#line 236 "zunmrz.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 237 "zunmrz.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 238 "zunmrz.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 242 "zunmrz.f"
    if (left) {
#line 243 "zunmrz.f"
	nq = *m;
#line 244 "zunmrz.f"
	nw = max(1,*n);
#line 245 "zunmrz.f"
    } else {
#line 246 "zunmrz.f"
	nq = *n;
#line 247 "zunmrz.f"
	nw = max(1,*m);
#line 248 "zunmrz.f"
    }
#line 249 "zunmrz.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 250 "zunmrz.f"
	*info = -1;
#line 251 "zunmrz.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 252 "zunmrz.f"
	*info = -2;
#line 253 "zunmrz.f"
    } else if (*m < 0) {
#line 254 "zunmrz.f"
	*info = -3;
#line 255 "zunmrz.f"
    } else if (*n < 0) {
#line 256 "zunmrz.f"
	*info = -4;
#line 257 "zunmrz.f"
    } else if (*k < 0 || *k > nq) {
#line 258 "zunmrz.f"
	*info = -5;
#line 259 "zunmrz.f"
    } else if (*l < 0 || left && *l > *m || ! left && *l > *n) {
#line 261 "zunmrz.f"
	*info = -6;
#line 262 "zunmrz.f"
    } else if (*lda < max(1,*k)) {
#line 263 "zunmrz.f"
	*info = -8;
#line 264 "zunmrz.f"
    } else if (*ldc < max(1,*m)) {
#line 265 "zunmrz.f"
	*info = -11;
#line 266 "zunmrz.f"
    }

#line 268 "zunmrz.f"
    if (*info == 0) {
#line 269 "zunmrz.f"
	if (*m == 0 || *n == 0) {
#line 270 "zunmrz.f"
	    lwkopt = 1;
#line 271 "zunmrz.f"
	} else {

/*           Determine the block size.  NB may be at most NBMAX, where */
/*           NBMAX is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 276 "zunmrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 276 "zunmrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 276 "zunmrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 276 "zunmrz.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 276 "zunmrz.f"
	    nb = min(i__1,i__2);
#line 278 "zunmrz.f"
	    lwkopt = nw * nb;
#line 279 "zunmrz.f"
	}
#line 280 "zunmrz.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 282 "zunmrz.f"
	if (*lwork < max(1,nw) && ! lquery) {
#line 283 "zunmrz.f"
	    *info = -13;
#line 284 "zunmrz.f"
	}
#line 285 "zunmrz.f"
    }

#line 287 "zunmrz.f"
    if (*info != 0) {
#line 288 "zunmrz.f"
	i__1 = -(*info);
#line 288 "zunmrz.f"
	xerbla_("ZUNMRZ", &i__1, (ftnlen)6);
#line 289 "zunmrz.f"
	return 0;
#line 290 "zunmrz.f"
    } else if (lquery) {
#line 291 "zunmrz.f"
	return 0;
#line 292 "zunmrz.f"
    }

/*     Quick return if possible */

#line 296 "zunmrz.f"
    if (*m == 0 || *n == 0) {
#line 297 "zunmrz.f"
	return 0;
#line 298 "zunmrz.f"
    }

/*     Determine the block size.  NB may be at most NBMAX, where NBMAX */
/*     is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 303 "zunmrz.f"
    i__3[0] = 1, a__1[0] = side;
#line 303 "zunmrz.f"
    i__3[1] = 1, a__1[1] = trans;
#line 303 "zunmrz.f"
    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 303 "zunmrz.f"
    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", ch__1, m, n, k, &c_n1, (ftnlen)
	    6, (ftnlen)2);
#line 303 "zunmrz.f"
    nb = min(i__1,i__2);
#line 305 "zunmrz.f"
    nbmin = 2;
#line 306 "zunmrz.f"
    ldwork = nw;
#line 307 "zunmrz.f"
    if (nb > 1 && nb < *k) {
#line 308 "zunmrz.f"
	iws = nw * nb;
#line 309 "zunmrz.f"
	if (*lwork < iws) {
#line 310 "zunmrz.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 311 "zunmrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 311 "zunmrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 311 "zunmrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 311 "zunmrz.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "ZUNMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 311 "zunmrz.f"
	    nbmin = max(i__1,i__2);
#line 313 "zunmrz.f"
	}
#line 314 "zunmrz.f"
    } else {
#line 315 "zunmrz.f"
	iws = nw;
#line 316 "zunmrz.f"
    }

#line 318 "zunmrz.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 322 "zunmrz.f"
	zunmr3_(side, trans, m, n, k, l, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 324 "zunmrz.f"
    } else {

/*        Use blocked code */

#line 328 "zunmrz.f"
	if (left && ! notran || ! left && notran) {
#line 330 "zunmrz.f"
	    i1 = 1;
#line 331 "zunmrz.f"
	    i2 = *k;
#line 332 "zunmrz.f"
	    i3 = nb;
#line 333 "zunmrz.f"
	} else {
#line 334 "zunmrz.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 335 "zunmrz.f"
	    i2 = 1;
#line 336 "zunmrz.f"
	    i3 = -nb;
#line 337 "zunmrz.f"
	}

#line 339 "zunmrz.f"
	if (left) {
#line 340 "zunmrz.f"
	    ni = *n;
#line 341 "zunmrz.f"
	    jc = 1;
#line 342 "zunmrz.f"
	    ja = *m - *l + 1;
#line 343 "zunmrz.f"
	} else {
#line 344 "zunmrz.f"
	    mi = *m;
#line 345 "zunmrz.f"
	    ic = 1;
#line 346 "zunmrz.f"
	    ja = *n - *l + 1;
#line 347 "zunmrz.f"
	}

#line 349 "zunmrz.f"
	if (notran) {
#line 350 "zunmrz.f"
	    *(unsigned char *)transt = 'C';
#line 351 "zunmrz.f"
	} else {
#line 352 "zunmrz.f"
	    *(unsigned char *)transt = 'N';
#line 353 "zunmrz.f"
	}

#line 355 "zunmrz.f"
	i__1 = i2;
#line 355 "zunmrz.f"
	i__2 = i3;
#line 355 "zunmrz.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 356 "zunmrz.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 356 "zunmrz.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 361 "zunmrz.f"
	    zlarzt_("Backward", "Rowwise", l, &ib, &a[i__ + ja * a_dim1], lda,
		     &tau[i__], t, &c__65, (ftnlen)8, (ftnlen)7);

#line 364 "zunmrz.f"
	    if (left) {

/*              H or H**H is applied to C(i:m,1:n) */

#line 368 "zunmrz.f"
		mi = *m - i__ + 1;
#line 369 "zunmrz.f"
		ic = i__;
#line 370 "zunmrz.f"
	    } else {

/*              H or H**H is applied to C(1:m,i:n) */

#line 374 "zunmrz.f"
		ni = *n - i__ + 1;
#line 375 "zunmrz.f"
		jc = i__;
#line 376 "zunmrz.f"
	    }

/*           Apply H or H**H */

#line 380 "zunmrz.f"
	    zlarzb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, l, &a[
		    i__ + ja * a_dim1], lda, t, &c__65, &c__[ic + jc * c_dim1]
		    , ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8,
		     (ftnlen)7);
#line 383 "zunmrz.f"
/* L10: */
#line 383 "zunmrz.f"
	}

#line 385 "zunmrz.f"
    }

#line 387 "zunmrz.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 389 "zunmrz.f"
    return 0;

/*     End of ZUNMRZ */

} /* zunmrz_ */

