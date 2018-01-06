#line 1 "cunmrz.f"
/* cunmrz.f -- translated by f2c (version 20100827).
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

#line 1 "cunmrz.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b CUNMRZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMRZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmrz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmrz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmrz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, L, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNMRZ overwrites the general complex M-by-N matrix C with */
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
/* > as returned by CTZRZF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CTZRZF in the last k rows of its array argument A. */
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
/* >          TAU is COMPLEX array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by CTZRZF. */
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
/* Subroutine */ int cunmrz_(char *side, char *trans, integer *m, integer *n, 
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
    extern /* Subroutine */ int cunmr3_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), clarzb_(char *, char *, char *, char *, integer *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), clarzt_(
	    char *, char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen);
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

#line 231 "cunmrz.f"
    /* Parameter adjustments */
#line 231 "cunmrz.f"
    a_dim1 = *lda;
#line 231 "cunmrz.f"
    a_offset = 1 + a_dim1;
#line 231 "cunmrz.f"
    a -= a_offset;
#line 231 "cunmrz.f"
    --tau;
#line 231 "cunmrz.f"
    c_dim1 = *ldc;
#line 231 "cunmrz.f"
    c_offset = 1 + c_dim1;
#line 231 "cunmrz.f"
    c__ -= c_offset;
#line 231 "cunmrz.f"
    --work;
#line 231 "cunmrz.f"

#line 231 "cunmrz.f"
    /* Function Body */
#line 231 "cunmrz.f"
    *info = 0;
#line 232 "cunmrz.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 233 "cunmrz.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 234 "cunmrz.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 238 "cunmrz.f"
    if (left) {
#line 239 "cunmrz.f"
	nq = *m;
#line 240 "cunmrz.f"
	nw = max(1,*n);
#line 241 "cunmrz.f"
    } else {
#line 242 "cunmrz.f"
	nq = *n;
#line 243 "cunmrz.f"
	nw = max(1,*m);
#line 244 "cunmrz.f"
    }
#line 245 "cunmrz.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 246 "cunmrz.f"
	*info = -1;
#line 247 "cunmrz.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 248 "cunmrz.f"
	*info = -2;
#line 249 "cunmrz.f"
    } else if (*m < 0) {
#line 250 "cunmrz.f"
	*info = -3;
#line 251 "cunmrz.f"
    } else if (*n < 0) {
#line 252 "cunmrz.f"
	*info = -4;
#line 253 "cunmrz.f"
    } else if (*k < 0 || *k > nq) {
#line 254 "cunmrz.f"
	*info = -5;
#line 255 "cunmrz.f"
    } else if (*l < 0 || left && *l > *m || ! left && *l > *n) {
#line 257 "cunmrz.f"
	*info = -6;
#line 258 "cunmrz.f"
    } else if (*lda < max(1,*k)) {
#line 259 "cunmrz.f"
	*info = -8;
#line 260 "cunmrz.f"
    } else if (*ldc < max(1,*m)) {
#line 261 "cunmrz.f"
	*info = -11;
#line 262 "cunmrz.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 263 "cunmrz.f"
	*info = -13;
#line 264 "cunmrz.f"
    }

#line 266 "cunmrz.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 270 "cunmrz.f"
	if (*m == 0 || *n == 0) {
#line 271 "cunmrz.f"
	    lwkopt = 1;
#line 272 "cunmrz.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 273 "cunmrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 273 "cunmrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 273 "cunmrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 273 "cunmrz.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "CUNMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 273 "cunmrz.f"
	    nb = min(i__1,i__2);
#line 275 "cunmrz.f"
	    lwkopt = nw * nb + 4160;
#line 276 "cunmrz.f"
	}
#line 277 "cunmrz.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 278 "cunmrz.f"
    }

#line 280 "cunmrz.f"
    if (*info != 0) {
#line 281 "cunmrz.f"
	i__1 = -(*info);
#line 281 "cunmrz.f"
	xerbla_("CUNMRZ", &i__1, (ftnlen)6);
#line 282 "cunmrz.f"
	return 0;
#line 283 "cunmrz.f"
    } else if (lquery) {
#line 284 "cunmrz.f"
	return 0;
#line 285 "cunmrz.f"
    }

/*     Quick return if possible */

#line 289 "cunmrz.f"
    if (*m == 0 || *n == 0) {
#line 290 "cunmrz.f"
	return 0;
#line 291 "cunmrz.f"
    }

/*     Determine the block size. */

/* Computing MIN */
/* Writing concatenation */
#line 295 "cunmrz.f"
    i__3[0] = 1, a__1[0] = side;
#line 295 "cunmrz.f"
    i__3[1] = 1, a__1[1] = trans;
#line 295 "cunmrz.f"
    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 295 "cunmrz.f"
    i__1 = 64, i__2 = ilaenv_(&c__1, "CUNMRQ", ch__1, m, n, k, &c_n1, (ftnlen)
	    6, (ftnlen)2);
#line 295 "cunmrz.f"
    nb = min(i__1,i__2);
#line 297 "cunmrz.f"
    nbmin = 2;
#line 298 "cunmrz.f"
    ldwork = nw;
#line 299 "cunmrz.f"
    if (nb > 1 && nb < *k) {
#line 300 "cunmrz.f"
	if (*lwork < nw * nb + 4160) {
#line 301 "cunmrz.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 302 "cunmrz.f"
	    i__3[0] = 1, a__1[0] = side;
#line 302 "cunmrz.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 302 "cunmrz.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 302 "cunmrz.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CUNMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 302 "cunmrz.f"
	    nbmin = max(i__1,i__2);
#line 304 "cunmrz.f"
	}
#line 305 "cunmrz.f"
    }

#line 307 "cunmrz.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 311 "cunmrz.f"
	cunmr3_(side, trans, m, n, k, l, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 313 "cunmrz.f"
    } else {

/*        Use blocked code */

#line 317 "cunmrz.f"
	iwt = nw * nb + 1;
#line 318 "cunmrz.f"
	if (left && ! notran || ! left && notran) {
#line 320 "cunmrz.f"
	    i1 = 1;
#line 321 "cunmrz.f"
	    i2 = *k;
#line 322 "cunmrz.f"
	    i3 = nb;
#line 323 "cunmrz.f"
	} else {
#line 324 "cunmrz.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 325 "cunmrz.f"
	    i2 = 1;
#line 326 "cunmrz.f"
	    i3 = -nb;
#line 327 "cunmrz.f"
	}

#line 329 "cunmrz.f"
	if (left) {
#line 330 "cunmrz.f"
	    ni = *n;
#line 331 "cunmrz.f"
	    jc = 1;
#line 332 "cunmrz.f"
	    ja = *m - *l + 1;
#line 333 "cunmrz.f"
	} else {
#line 334 "cunmrz.f"
	    mi = *m;
#line 335 "cunmrz.f"
	    ic = 1;
#line 336 "cunmrz.f"
	    ja = *n - *l + 1;
#line 337 "cunmrz.f"
	}

#line 339 "cunmrz.f"
	if (notran) {
#line 340 "cunmrz.f"
	    *(unsigned char *)transt = 'C';
#line 341 "cunmrz.f"
	} else {
#line 342 "cunmrz.f"
	    *(unsigned char *)transt = 'N';
#line 343 "cunmrz.f"
	}

#line 345 "cunmrz.f"
	i__1 = i2;
#line 345 "cunmrz.f"
	i__2 = i3;
#line 345 "cunmrz.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 346 "cunmrz.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 346 "cunmrz.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 351 "cunmrz.f"
	    clarzt_("Backward", "Rowwise", l, &ib, &a[i__ + ja * a_dim1], lda,
		     &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)7);

#line 354 "cunmrz.f"
	    if (left) {

/*              H or H**H is applied to C(i:m,1:n) */

#line 358 "cunmrz.f"
		mi = *m - i__ + 1;
#line 359 "cunmrz.f"
		ic = i__;
#line 360 "cunmrz.f"
	    } else {

/*              H or H**H is applied to C(1:m,i:n) */

#line 364 "cunmrz.f"
		ni = *n - i__ + 1;
#line 365 "cunmrz.f"
		jc = i__;
#line 366 "cunmrz.f"
	    }

/*           Apply H or H**H */

#line 370 "cunmrz.f"
	    clarzb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, l, &a[
		    i__ + ja * a_dim1], lda, &work[iwt], &c__65, &c__[ic + jc 
		    * c_dim1], ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (
		    ftnlen)8, (ftnlen)7);
#line 373 "cunmrz.f"
/* L10: */
#line 373 "cunmrz.f"
	}

#line 375 "cunmrz.f"
    }

#line 377 "cunmrz.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 379 "cunmrz.f"
    return 0;

/*     End of CUNMRZ */

} /* cunmrz_ */

