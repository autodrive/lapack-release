#line 1 "cunmlq.f"
/* cunmlq.f -- translated by f2c (version 20100827).
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

#line 1 "cunmlq.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b CUNMLQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMLQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmlq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmlq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmlq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > CUNMLQ overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(k)**H . . . H(2)**H H(1)**H */
/* > */
/* > as returned by CGELQF. Q is of order M if SIDE = 'L' and of order N */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension */
/* >                               (LDA,M) if SIDE = 'L', */
/* >                               (LDA,N) if SIDE = 'R' */
/* >          The i-th row must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGELQF in the first k rows of its array argument A. */
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
/* >          reflector H(i), as returned by CGELQF. */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunmlq_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i__, i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iwt;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int cunml2_(char *, char *, integer *, integer *, 
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

#line 213 "cunmlq.f"
    /* Parameter adjustments */
#line 213 "cunmlq.f"
    a_dim1 = *lda;
#line 213 "cunmlq.f"
    a_offset = 1 + a_dim1;
#line 213 "cunmlq.f"
    a -= a_offset;
#line 213 "cunmlq.f"
    --tau;
#line 213 "cunmlq.f"
    c_dim1 = *ldc;
#line 213 "cunmlq.f"
    c_offset = 1 + c_dim1;
#line 213 "cunmlq.f"
    c__ -= c_offset;
#line 213 "cunmlq.f"
    --work;
#line 213 "cunmlq.f"

#line 213 "cunmlq.f"
    /* Function Body */
#line 213 "cunmlq.f"
    *info = 0;
#line 214 "cunmlq.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 215 "cunmlq.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 216 "cunmlq.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 220 "cunmlq.f"
    if (left) {
#line 221 "cunmlq.f"
	nq = *m;
#line 222 "cunmlq.f"
	nw = *n;
#line 223 "cunmlq.f"
    } else {
#line 224 "cunmlq.f"
	nq = *n;
#line 225 "cunmlq.f"
	nw = *m;
#line 226 "cunmlq.f"
    }
#line 227 "cunmlq.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 228 "cunmlq.f"
	*info = -1;
#line 229 "cunmlq.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 230 "cunmlq.f"
	*info = -2;
#line 231 "cunmlq.f"
    } else if (*m < 0) {
#line 232 "cunmlq.f"
	*info = -3;
#line 233 "cunmlq.f"
    } else if (*n < 0) {
#line 234 "cunmlq.f"
	*info = -4;
#line 235 "cunmlq.f"
    } else if (*k < 0 || *k > nq) {
#line 236 "cunmlq.f"
	*info = -5;
#line 237 "cunmlq.f"
    } else if (*lda < max(1,*k)) {
#line 238 "cunmlq.f"
	*info = -7;
#line 239 "cunmlq.f"
    } else if (*ldc < max(1,*m)) {
#line 240 "cunmlq.f"
	*info = -10;
#line 241 "cunmlq.f"
    } else if (*lwork < max(1,nw) && ! lquery) {
#line 242 "cunmlq.f"
	*info = -12;
#line 243 "cunmlq.f"
    }

#line 245 "cunmlq.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 249 "cunmlq.f"
	if (*m == 0 || *n == 0 || *k == 0) {
#line 250 "cunmlq.f"
	    lwkopt = 1;
#line 251 "cunmlq.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 252 "cunmlq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 252 "cunmlq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 252 "cunmlq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 252 "cunmlq.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "CUNMLQ", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 252 "cunmlq.f"
	    nb = min(i__1,i__2);
#line 254 "cunmlq.f"
	    lwkopt = max(1,nw) * nb + 4160;
#line 255 "cunmlq.f"
	}
#line 256 "cunmlq.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 257 "cunmlq.f"
    }

#line 259 "cunmlq.f"
    if (*info != 0) {
#line 260 "cunmlq.f"
	i__1 = -(*info);
#line 260 "cunmlq.f"
	xerbla_("CUNMLQ", &i__1, (ftnlen)6);
#line 261 "cunmlq.f"
	return 0;
#line 262 "cunmlq.f"
    } else if (lquery) {
#line 263 "cunmlq.f"
	return 0;
#line 264 "cunmlq.f"
    }

/*     Quick return if possible */

#line 268 "cunmlq.f"
    if (*m == 0 || *n == 0 || *k == 0) {
#line 269 "cunmlq.f"
	return 0;
#line 270 "cunmlq.f"
    }

/*     Determine the block size */

#line 274 "cunmlq.f"
    nbmin = 2;
#line 275 "cunmlq.f"
    ldwork = nw;
#line 276 "cunmlq.f"
    if (nb > 1 && nb < *k) {
#line 277 "cunmlq.f"
	if (*lwork < nw * nb + 4160) {
#line 278 "cunmlq.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 279 "cunmlq.f"
	    i__3[0] = 1, a__1[0] = side;
#line 279 "cunmlq.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 279 "cunmlq.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 279 "cunmlq.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CUNMLQ", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 279 "cunmlq.f"
	    nbmin = max(i__1,i__2);
#line 281 "cunmlq.f"
	}
#line 282 "cunmlq.f"
    }

#line 284 "cunmlq.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 288 "cunmlq.f"
	cunml2_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 290 "cunmlq.f"
    } else {

/*        Use blocked code */

#line 294 "cunmlq.f"
	iwt = nw * nb + 1;
#line 295 "cunmlq.f"
	if (left && notran || ! left && ! notran) {
#line 297 "cunmlq.f"
	    i1 = 1;
#line 298 "cunmlq.f"
	    i2 = *k;
#line 299 "cunmlq.f"
	    i3 = nb;
#line 300 "cunmlq.f"
	} else {
#line 301 "cunmlq.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 302 "cunmlq.f"
	    i2 = 1;
#line 303 "cunmlq.f"
	    i3 = -nb;
#line 304 "cunmlq.f"
	}

#line 306 "cunmlq.f"
	if (left) {
#line 307 "cunmlq.f"
	    ni = *n;
#line 308 "cunmlq.f"
	    jc = 1;
#line 309 "cunmlq.f"
	} else {
#line 310 "cunmlq.f"
	    mi = *m;
#line 311 "cunmlq.f"
	    ic = 1;
#line 312 "cunmlq.f"
	}

#line 314 "cunmlq.f"
	if (notran) {
#line 315 "cunmlq.f"
	    *(unsigned char *)transt = 'C';
#line 316 "cunmlq.f"
	} else {
#line 317 "cunmlq.f"
	    *(unsigned char *)transt = 'N';
#line 318 "cunmlq.f"
	}

#line 320 "cunmlq.f"
	i__1 = i2;
#line 320 "cunmlq.f"
	i__2 = i3;
#line 320 "cunmlq.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 321 "cunmlq.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 321 "cunmlq.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i) H(i+1) . . . H(i+ib-1) */

#line 326 "cunmlq.f"
	    i__4 = nq - i__ + 1;
#line 326 "cunmlq.f"
	    clarft_("Forward", "Rowwise", &i__4, &ib, &a[i__ + i__ * a_dim1], 
		    lda, &tau[i__], &work[iwt], &c__65, (ftnlen)7, (ftnlen)7);
#line 328 "cunmlq.f"
	    if (left) {

/*              H or H**H is applied to C(i:m,1:n) */

#line 332 "cunmlq.f"
		mi = *m - i__ + 1;
#line 333 "cunmlq.f"
		ic = i__;
#line 334 "cunmlq.f"
	    } else {

/*              H or H**H is applied to C(1:m,i:n) */

#line 338 "cunmlq.f"
		ni = *n - i__ + 1;
#line 339 "cunmlq.f"
		jc = i__;
#line 340 "cunmlq.f"
	    }

/*           Apply H or H**H */

#line 344 "cunmlq.f"
	    clarfb_(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &a[i__ 
		    + i__ * a_dim1], lda, &work[iwt], &c__65, &c__[ic + jc * 
		    c_dim1], ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (
		    ftnlen)7, (ftnlen)7);
#line 347 "cunmlq.f"
/* L10: */
#line 347 "cunmlq.f"
	}
#line 348 "cunmlq.f"
    }
#line 349 "cunmlq.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 350 "cunmlq.f"
    return 0;

/*     End of CUNMLQ */

} /* cunmlq_ */

