#line 1 "dormql.f"
/* dormql.f -- translated by f2c (version 20100827).
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

#line 1 "dormql.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b DORMQL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORMQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormql.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormql.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormql.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORMQL overwrites the general real M-by-N matrix C with */
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
/* > as returned by DGEQLF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DGEQLF in the last k columns of its array argument A. */
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
/* >          TAU is DOUBLE PRECISION array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGEQLF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dormql_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i__, i1, i2, i3, ib, nb, mi, ni, nq, nw, iwt;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int dorm2l_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen), dlarfb_(char 
	    *, char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), dlarft_(char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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

#line 210 "dormql.f"
    /* Parameter adjustments */
#line 210 "dormql.f"
    a_dim1 = *lda;
#line 210 "dormql.f"
    a_offset = 1 + a_dim1;
#line 210 "dormql.f"
    a -= a_offset;
#line 210 "dormql.f"
    --tau;
#line 210 "dormql.f"
    c_dim1 = *ldc;
#line 210 "dormql.f"
    c_offset = 1 + c_dim1;
#line 210 "dormql.f"
    c__ -= c_offset;
#line 210 "dormql.f"
    --work;
#line 210 "dormql.f"

#line 210 "dormql.f"
    /* Function Body */
#line 210 "dormql.f"
    *info = 0;
#line 211 "dormql.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 212 "dormql.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 213 "dormql.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 217 "dormql.f"
    if (left) {
#line 218 "dormql.f"
	nq = *m;
#line 219 "dormql.f"
	nw = max(1,*n);
#line 220 "dormql.f"
    } else {
#line 221 "dormql.f"
	nq = *n;
#line 222 "dormql.f"
	nw = max(1,*m);
#line 223 "dormql.f"
    }
#line 224 "dormql.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 225 "dormql.f"
	*info = -1;
#line 226 "dormql.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 227 "dormql.f"
	*info = -2;
#line 228 "dormql.f"
    } else if (*m < 0) {
#line 229 "dormql.f"
	*info = -3;
#line 230 "dormql.f"
    } else if (*n < 0) {
#line 231 "dormql.f"
	*info = -4;
#line 232 "dormql.f"
    } else if (*k < 0 || *k > nq) {
#line 233 "dormql.f"
	*info = -5;
#line 234 "dormql.f"
    } else if (*lda < max(1,nq)) {
#line 235 "dormql.f"
	*info = -7;
#line 236 "dormql.f"
    } else if (*ldc < max(1,*m)) {
#line 237 "dormql.f"
	*info = -10;
#line 238 "dormql.f"
    } else if (*lwork < nw && ! lquery) {
#line 239 "dormql.f"
	*info = -12;
#line 240 "dormql.f"
    }

#line 242 "dormql.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 246 "dormql.f"
	if (*m == 0 || *n == 0) {
#line 247 "dormql.f"
	    lwkopt = 1;
#line 248 "dormql.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 249 "dormql.f"
	    i__3[0] = 1, a__1[0] = side;
#line 249 "dormql.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 249 "dormql.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 249 "dormql.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQL", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 249 "dormql.f"
	    nb = min(i__1,i__2);
#line 251 "dormql.f"
	    lwkopt = nw * nb + 4160;
#line 252 "dormql.f"
	}
#line 253 "dormql.f"
	work[1] = (doublereal) lwkopt;
#line 254 "dormql.f"
    }

#line 256 "dormql.f"
    if (*info != 0) {
#line 257 "dormql.f"
	i__1 = -(*info);
#line 257 "dormql.f"
	xerbla_("DORMQL", &i__1, (ftnlen)6);
#line 258 "dormql.f"
	return 0;
#line 259 "dormql.f"
    } else if (lquery) {
#line 260 "dormql.f"
	return 0;
#line 261 "dormql.f"
    }

/*     Quick return if possible */

#line 265 "dormql.f"
    if (*m == 0 || *n == 0) {
#line 266 "dormql.f"
	return 0;
#line 267 "dormql.f"
    }

#line 269 "dormql.f"
    nbmin = 2;
#line 270 "dormql.f"
    ldwork = nw;
#line 271 "dormql.f"
    if (nb > 1 && nb < *k) {
#line 272 "dormql.f"
	if (*lwork < nw * nb + 4160) {
#line 273 "dormql.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 274 "dormql.f"
	    i__3[0] = 1, a__1[0] = side;
#line 274 "dormql.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 274 "dormql.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 274 "dormql.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMQL", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 274 "dormql.f"
	    nbmin = max(i__1,i__2);
#line 276 "dormql.f"
	}
#line 277 "dormql.f"
    }

#line 279 "dormql.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 283 "dormql.f"
	dorm2l_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 285 "dormql.f"
    } else {

/*        Use blocked code */

#line 289 "dormql.f"
	iwt = nw * nb + 1;
#line 290 "dormql.f"
	if (left && notran || ! left && ! notran) {
#line 292 "dormql.f"
	    i1 = 1;
#line 293 "dormql.f"
	    i2 = *k;
#line 294 "dormql.f"
	    i3 = nb;
#line 295 "dormql.f"
	} else {
#line 296 "dormql.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 297 "dormql.f"
	    i2 = 1;
#line 298 "dormql.f"
	    i3 = -nb;
#line 299 "dormql.f"
	}

#line 301 "dormql.f"
	if (left) {
#line 302 "dormql.f"
	    ni = *n;
#line 303 "dormql.f"
	} else {
#line 304 "dormql.f"
	    mi = *m;
#line 305 "dormql.f"
	}

#line 307 "dormql.f"
	i__1 = i2;
#line 307 "dormql.f"
	i__2 = i3;
#line 307 "dormql.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 308 "dormql.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 308 "dormql.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 313 "dormql.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 313 "dormql.f"
	    dlarft_("Backward", "Columnwise", &i__4, &ib, &a[i__ * a_dim1 + 1]
		    , lda, &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)
		    10);
#line 315 "dormql.f"
	    if (left) {

/*              H or H**T is applied to C(1:m-k+i+ib-1,1:n) */

#line 319 "dormql.f"
		mi = *m - *k + i__ + ib - 1;
#line 320 "dormql.f"
	    } else {

/*              H or H**T is applied to C(1:m,1:n-k+i+ib-1) */

#line 324 "dormql.f"
		ni = *n - *k + i__ + ib - 1;
#line 325 "dormql.f"
	    }

/*           Apply H or H**T */

#line 329 "dormql.f"
	    dlarfb_(side, trans, "Backward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ * a_dim1 + 1], lda, &work[iwt], &c__65, &c__[c_offset]
		    , ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8,
		     (ftnlen)10);
#line 332 "dormql.f"
/* L10: */
#line 332 "dormql.f"
	}
#line 333 "dormql.f"
    }
#line 334 "dormql.f"
    work[1] = (doublereal) lwkopt;
#line 335 "dormql.f"
    return 0;

/*     End of DORMQL */

} /* dormql_ */

