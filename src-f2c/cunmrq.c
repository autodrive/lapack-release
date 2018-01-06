#line 1 "cunmrq.f"
/* cunmrq.f -- translated by f2c (version 20100827).
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

#line 1 "cunmrq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b CUNMRQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMRQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmrq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmrq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmrq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNMRQ overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(1)**H H(2)**H . . . H(k)**H */
/* > */
/* > as returned by CGERQF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          = 'C':  Transpose, apply Q**H. */
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
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGERQF in the last k rows of its array argument A. */
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
/* >          reflector H(i), as returned by CGERQF. */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunmrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
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
    static doublecomplex t[4160]	/* was [65][64] */;
    static integer i1, i2, i3, ib, nb, mi, ni, nq, nw, iws;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int cunmr2_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), clarfb_(char *, char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublecomplex *, integer *
	    , ftnlen, ftnlen, ftnlen, ftnlen), clarft_(char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, ftnlen, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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

#line 217 "cunmrq.f"
    /* Parameter adjustments */
#line 217 "cunmrq.f"
    a_dim1 = *lda;
#line 217 "cunmrq.f"
    a_offset = 1 + a_dim1;
#line 217 "cunmrq.f"
    a -= a_offset;
#line 217 "cunmrq.f"
    --tau;
#line 217 "cunmrq.f"
    c_dim1 = *ldc;
#line 217 "cunmrq.f"
    c_offset = 1 + c_dim1;
#line 217 "cunmrq.f"
    c__ -= c_offset;
#line 217 "cunmrq.f"
    --work;
#line 217 "cunmrq.f"

#line 217 "cunmrq.f"
    /* Function Body */
#line 217 "cunmrq.f"
    *info = 0;
#line 218 "cunmrq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 219 "cunmrq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 220 "cunmrq.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 224 "cunmrq.f"
    if (left) {
#line 225 "cunmrq.f"
	nq = *m;
#line 226 "cunmrq.f"
	nw = max(1,*n);
#line 227 "cunmrq.f"
    } else {
#line 228 "cunmrq.f"
	nq = *n;
#line 229 "cunmrq.f"
	nw = max(1,*m);
#line 230 "cunmrq.f"
    }
#line 231 "cunmrq.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 232 "cunmrq.f"
	*info = -1;
#line 233 "cunmrq.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 234 "cunmrq.f"
	*info = -2;
#line 235 "cunmrq.f"
    } else if (*m < 0) {
#line 236 "cunmrq.f"
	*info = -3;
#line 237 "cunmrq.f"
    } else if (*n < 0) {
#line 238 "cunmrq.f"
	*info = -4;
#line 239 "cunmrq.f"
    } else if (*k < 0 || *k > nq) {
#line 240 "cunmrq.f"
	*info = -5;
#line 241 "cunmrq.f"
    } else if (*lda < max(1,*k)) {
#line 242 "cunmrq.f"
	*info = -7;
#line 243 "cunmrq.f"
    } else if (*ldc < max(1,*m)) {
#line 244 "cunmrq.f"
	*info = -10;
#line 245 "cunmrq.f"
    }

#line 247 "cunmrq.f"
    if (*info == 0) {
#line 248 "cunmrq.f"
	if (*m == 0 || *n == 0) {
#line 249 "cunmrq.f"
	    lwkopt = 1;
#line 250 "cunmrq.f"
	} else {

/*           Determine the block size.  NB may be at most NBMAX, where */
/*           NBMAX is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
#line 255 "cunmrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 255 "cunmrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 255 "cunmrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 255 "cunmrq.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "CUNMRQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 255 "cunmrq.f"
	    nb = min(i__1,i__2);
#line 257 "cunmrq.f"
	    lwkopt = nw * nb;
#line 258 "cunmrq.f"
	}
#line 259 "cunmrq.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 261 "cunmrq.f"
	if (*lwork < nw && ! lquery) {
#line 262 "cunmrq.f"
	    *info = -12;
#line 263 "cunmrq.f"
	}
#line 264 "cunmrq.f"
    }

#line 266 "cunmrq.f"
    if (*info != 0) {
#line 267 "cunmrq.f"
	i__1 = -(*info);
#line 267 "cunmrq.f"
	xerbla_("CUNMRQ", &i__1, (ftnlen)6);
#line 268 "cunmrq.f"
	return 0;
#line 269 "cunmrq.f"
    } else if (lquery) {
#line 270 "cunmrq.f"
	return 0;
#line 271 "cunmrq.f"
    }

/*     Quick return if possible */

#line 275 "cunmrq.f"
    if (*m == 0 || *n == 0) {
#line 276 "cunmrq.f"
	return 0;
#line 277 "cunmrq.f"
    }

#line 279 "cunmrq.f"
    nbmin = 2;
#line 280 "cunmrq.f"
    ldwork = nw;
#line 281 "cunmrq.f"
    if (nb > 1 && nb < *k) {
#line 282 "cunmrq.f"
	iws = nw * nb;
#line 283 "cunmrq.f"
	if (*lwork < iws) {
#line 284 "cunmrq.f"
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 285 "cunmrq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 285 "cunmrq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 285 "cunmrq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 285 "cunmrq.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CUNMRQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 285 "cunmrq.f"
	    nbmin = max(i__1,i__2);
#line 287 "cunmrq.f"
	}
#line 288 "cunmrq.f"
    } else {
#line 289 "cunmrq.f"
	iws = nw;
#line 290 "cunmrq.f"
    }

#line 292 "cunmrq.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 296 "cunmrq.f"
	cunmr2_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 298 "cunmrq.f"
    } else {

/*        Use blocked code */

#line 302 "cunmrq.f"
	if (left && ! notran || ! left && notran) {
#line 304 "cunmrq.f"
	    i1 = 1;
#line 305 "cunmrq.f"
	    i2 = *k;
#line 306 "cunmrq.f"
	    i3 = nb;
#line 307 "cunmrq.f"
	} else {
#line 308 "cunmrq.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 309 "cunmrq.f"
	    i2 = 1;
#line 310 "cunmrq.f"
	    i3 = -nb;
#line 311 "cunmrq.f"
	}

#line 313 "cunmrq.f"
	if (left) {
#line 314 "cunmrq.f"
	    ni = *n;
#line 315 "cunmrq.f"
	} else {
#line 316 "cunmrq.f"
	    mi = *m;
#line 317 "cunmrq.f"
	}

#line 319 "cunmrq.f"
	if (notran) {
#line 320 "cunmrq.f"
	    *(unsigned char *)transt = 'C';
#line 321 "cunmrq.f"
	} else {
#line 322 "cunmrq.f"
	    *(unsigned char *)transt = 'N';
#line 323 "cunmrq.f"
	}

#line 325 "cunmrq.f"
	i__1 = i2;
#line 325 "cunmrq.f"
	i__2 = i3;
#line 325 "cunmrq.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 326 "cunmrq.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 326 "cunmrq.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 331 "cunmrq.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 331 "cunmrq.f"
	    clarft_("Backward", "Rowwise", &i__4, &ib, &a[i__ + a_dim1], lda, 
		    &tau[i__], t, &c__65, (ftnlen)8, (ftnlen)7);
#line 333 "cunmrq.f"
	    if (left) {

/*              H or H**H is applied to C(1:m-k+i+ib-1,1:n) */

#line 337 "cunmrq.f"
		mi = *m - *k + i__ + ib - 1;
#line 338 "cunmrq.f"
	    } else {

/*              H or H**H is applied to C(1:m,1:n-k+i+ib-1) */

#line 342 "cunmrq.f"
		ni = *n - *k + i__ + ib - 1;
#line 343 "cunmrq.f"
	    }

/*           Apply H or H**H */

#line 347 "cunmrq.f"
	    clarfb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, &a[
		    i__ + a_dim1], lda, t, &c__65, &c__[c_offset], ldc, &work[
		    1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8, (ftnlen)7);
#line 350 "cunmrq.f"
/* L10: */
#line 350 "cunmrq.f"
	}
#line 351 "cunmrq.f"
    }
#line 352 "cunmrq.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 353 "cunmrq.f"
    return 0;

/*     End of CUNMRQ */

} /* cunmrq_ */

