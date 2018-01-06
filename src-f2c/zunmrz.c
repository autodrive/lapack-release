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
    static integer i__, i1, i2, i3, ib, ic, ja, jc, nb, mi, ni, nq, nw, iwt;
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

#line 231 "zunmrz.f"
    /* Parameter adjustments */
#line 231 "zunmrz.f"
    a_dim1 = *lda;
#line 231 "zunmrz.f"
    a_offset = 1 + a_dim1;
#line 231 "zunmrz.f"
    a -= a_offset;
#line 231 "zunmrz.f"
    --tau;
#line 231 "zunmrz.f"
    c_dim1 = *ldc;
#line 231 "zunmrz.f"
    c_offset = 1 + c_dim1;
#line 231 "zunmrz.f"
    c__ -= c_offset;
#line 231 "zunmrz.f"
    --work;
#line 231 "zunmrz.f"

#line 231 "zunmrz.f"
    /* Function Body */
#line 231 "zunmrz.f"
    *info = 0;
#line 232 "zunmrz.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 233 "zunmrz.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 234 "zunmrz.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 238 "zunmrz.f"
    if (left) {
#line 239 "zunmrz.f"
	nq = *m;
#line 240 "zunmrz.f"
	nw = max(1,*n);
#line 241 "zunmrz.f"
    } else {
#line 242 "zunmrz.f"
	nq = *n;
#line 243 "zunmrz.f"
	nw = max(1,*m);
#line 244 "zunmrz.f"
    }
#line 245 "zunmrz.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 246 "zunmrz.f"
	*info = -1;
#line 247 "zunmrz.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 248 "zunmrz.f"
	*info = -2;
#line 249 "zunmrz.f"
    } else if (*m < 0) {
#line 250 "zunmrz.f"
	*info = -3;
#line 251 "zunmrz.f"
    } else if (*n < 0) {
#line 252 "zunmrz.f"
	*info = -4;
#line 253 "zunmrz.f"
    } else if (*k < 0 || *k > nq) {
#line 254 "zunmrz.f"
	*info = -5;
#line 255 "zunmrz.f"
    } else if (*l < 0 || left && *l > *m || ! left && *l > *n) {
#line 257 "zunmrz.f"
	*info = -6;
#line 258 "zunmrz.f"
    } else if (*lda < max(1,*k)) {
#line 259 "zunmrz.f"
	*info = -8;
#line 260 "zunmrz.f"
    } else if (*ldc < max(1,*m)) {
#line 261 "zunmrz.f"
	*info = -11;
#line 262 "zunmrz.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 263 "zunmrz.f"
	*info = -13;
#line 264 "zunmrz.f"
    }

#line 266 "zunmrz.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 270 "zunmrz.f"
	if (*m == 0 || *n == 0) {
#line 271 "zunmrz.f"
	    lwkopt = 1;
#line 272 "zunmrz.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 273 "zunmrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 273 "zunmrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 273 "zunmrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 273 "zunmrz.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 273 "zunmrz.f"
	    nb = min(i__1,i__2);
#line 275 "zunmrz.f"
	    lwkopt = nw * nb + 4160;
#line 276 "zunmrz.f"
	}
#line 277 "zunmrz.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 278 "zunmrz.f"
    }

#line 280 "zunmrz.f"
    if (*info != 0) {
#line 281 "zunmrz.f"
	i__1 = -(*info);
#line 281 "zunmrz.f"
	xerbla_("ZUNMRZ", &i__1, (ftnlen)6);
#line 282 "zunmrz.f"
	return 0;
#line 283 "zunmrz.f"
    } else if (lquery) {
#line 284 "zunmrz.f"
	return 0;
#line 285 "zunmrz.f"
    }

/*     Quick return if possible */

#line 289 "zunmrz.f"
    if (*m == 0 || *n == 0) {
#line 290 "zunmrz.f"
	return 0;
#line 291 "zunmrz.f"
    }

/*     Determine the block size.  NB may be at most NBMAX, where NBMAX */
/*     is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 296 "zunmrz.f"
    i__3[0] = 1, a__1[0] = side;
#line 296 "zunmrz.f"
    i__3[1] = 1, a__1[1] = trans;
#line 296 "zunmrz.f"
    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 296 "zunmrz.f"
    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", ch__1, m, n, k, &c_n1, (ftnlen)
	    6, (ftnlen)2);
#line 296 "zunmrz.f"
    nb = min(i__1,i__2);
#line 298 "zunmrz.f"
    nbmin = 2;
#line 299 "zunmrz.f"
    ldwork = nw;
#line 300 "zunmrz.f"
    if (nb > 1 && nb < *k) {
#line 301 "zunmrz.f"
	if (*lwork < nw * nb + 4160) {
#line 302 "zunmrz.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 303 "zunmrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 303 "zunmrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 303 "zunmrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 303 "zunmrz.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "ZUNMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 303 "zunmrz.f"
	    nbmin = max(i__1,i__2);
#line 305 "zunmrz.f"
	}
#line 306 "zunmrz.f"
    }

#line 308 "zunmrz.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 312 "zunmrz.f"
	zunmr3_(side, trans, m, n, k, l, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 314 "zunmrz.f"
    } else {

/*        Use blocked code */

#line 318 "zunmrz.f"
	iwt = nw * nb + 1;
#line 319 "zunmrz.f"
	if (left && ! notran || ! left && notran) {
#line 321 "zunmrz.f"
	    i1 = 1;
#line 322 "zunmrz.f"
	    i2 = *k;
#line 323 "zunmrz.f"
	    i3 = nb;
#line 324 "zunmrz.f"
	} else {
#line 325 "zunmrz.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 326 "zunmrz.f"
	    i2 = 1;
#line 327 "zunmrz.f"
	    i3 = -nb;
#line 328 "zunmrz.f"
	}

#line 330 "zunmrz.f"
	if (left) {
#line 331 "zunmrz.f"
	    ni = *n;
#line 332 "zunmrz.f"
	    jc = 1;
#line 333 "zunmrz.f"
	    ja = *m - *l + 1;
#line 334 "zunmrz.f"
	} else {
#line 335 "zunmrz.f"
	    mi = *m;
#line 336 "zunmrz.f"
	    ic = 1;
#line 337 "zunmrz.f"
	    ja = *n - *l + 1;
#line 338 "zunmrz.f"
	}

#line 340 "zunmrz.f"
	if (notran) {
#line 341 "zunmrz.f"
	    *(unsigned char *)transt = 'C';
#line 342 "zunmrz.f"
	} else {
#line 343 "zunmrz.f"
	    *(unsigned char *)transt = 'N';
#line 344 "zunmrz.f"
	}

#line 346 "zunmrz.f"
	i__1 = i2;
#line 346 "zunmrz.f"
	i__2 = i3;
#line 346 "zunmrz.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 347 "zunmrz.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 347 "zunmrz.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 352 "zunmrz.f"
	    zlarzt_("Backward", "Rowwise", l, &ib, &a[i__ + ja * a_dim1], lda,
		     &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)7);

#line 355 "zunmrz.f"
	    if (left) {

/*              H or H**H is applied to C(i:m,1:n) */

#line 359 "zunmrz.f"
		mi = *m - i__ + 1;
#line 360 "zunmrz.f"
		ic = i__;
#line 361 "zunmrz.f"
	    } else {

/*              H or H**H is applied to C(1:m,i:n) */

#line 365 "zunmrz.f"
		ni = *n - i__ + 1;
#line 366 "zunmrz.f"
		jc = i__;
#line 367 "zunmrz.f"
	    }

/*           Apply H or H**H */

#line 371 "zunmrz.f"
	    zlarzb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, l, &a[
		    i__ + ja * a_dim1], lda, &work[iwt], &c__65, &c__[ic + jc 
		    * c_dim1], ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (
		    ftnlen)8, (ftnlen)7);
#line 374 "zunmrz.f"
/* L10: */
#line 374 "zunmrz.f"
	}

#line 376 "zunmrz.f"
    }

#line 378 "zunmrz.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 380 "zunmrz.f"
    return 0;

/*     End of ZUNMRZ */

} /* zunmrz_ */

