#line 1 "cunmql.f"
/* cunmql.f -- translated by f2c (version 20100827).
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

#line 1 "cunmql.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* > \brief \b CUNMQL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmql.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmql.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmql.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > CUNMQL overwrites the general complex M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by CGEQLF. Q is of order M if SIDE = 'L' and of order N */
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
/* >          A is COMPLEX array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          CGEQLF in the last k columns of its array argument A. */
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
/* >          TAU is COMPLEX array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by CGEQLF. */
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

/*  ===================================================================== */
/* Subroutine */ int cunmql_(char *side, char *trans, integer *m, integer *n, 
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
    static integer i__, i1, i2, i3, ib, nb, mi, ni, nq, nw, iwt;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    extern /* Subroutine */ int cunm2l_(char *, char *, integer *, integer *, 
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
    static integer ldwork, lwkopt;
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

#line 212 "cunmql.f"
    /* Parameter adjustments */
#line 212 "cunmql.f"
    a_dim1 = *lda;
#line 212 "cunmql.f"
    a_offset = 1 + a_dim1;
#line 212 "cunmql.f"
    a -= a_offset;
#line 212 "cunmql.f"
    --tau;
#line 212 "cunmql.f"
    c_dim1 = *ldc;
#line 212 "cunmql.f"
    c_offset = 1 + c_dim1;
#line 212 "cunmql.f"
    c__ -= c_offset;
#line 212 "cunmql.f"
    --work;
#line 212 "cunmql.f"

#line 212 "cunmql.f"
    /* Function Body */
#line 212 "cunmql.f"
    *info = 0;
#line 213 "cunmql.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 214 "cunmql.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 215 "cunmql.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

#line 219 "cunmql.f"
    if (left) {
#line 220 "cunmql.f"
	nq = *m;
#line 221 "cunmql.f"
	nw = max(1,*n);
#line 222 "cunmql.f"
    } else {
#line 223 "cunmql.f"
	nq = *n;
#line 224 "cunmql.f"
	nw = max(1,*m);
#line 225 "cunmql.f"
    }
#line 226 "cunmql.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 227 "cunmql.f"
	*info = -1;
#line 228 "cunmql.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 229 "cunmql.f"
	*info = -2;
#line 230 "cunmql.f"
    } else if (*m < 0) {
#line 231 "cunmql.f"
	*info = -3;
#line 232 "cunmql.f"
    } else if (*n < 0) {
#line 233 "cunmql.f"
	*info = -4;
#line 234 "cunmql.f"
    } else if (*k < 0 || *k > nq) {
#line 235 "cunmql.f"
	*info = -5;
#line 236 "cunmql.f"
    } else if (*lda < max(1,nq)) {
#line 237 "cunmql.f"
	*info = -7;
#line 238 "cunmql.f"
    } else if (*ldc < max(1,*m)) {
#line 239 "cunmql.f"
	*info = -10;
#line 240 "cunmql.f"
    } else if (*lwork < nw && ! lquery) {
#line 241 "cunmql.f"
	*info = -12;
#line 242 "cunmql.f"
    }

#line 244 "cunmql.f"
    if (*info == 0) {

/*        Compute the workspace requirements */

#line 248 "cunmql.f"
	if (*m == 0 || *n == 0) {
#line 249 "cunmql.f"
	    lwkopt = 1;
#line 250 "cunmql.f"
	} else {
/* Computing MIN */
/* Writing concatenation */
#line 251 "cunmql.f"
	    i__3[0] = 1, a__1[0] = side;
#line 251 "cunmql.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 251 "cunmql.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 251 "cunmql.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "CUNMQL", ch__1, m, n, k, &c_n1, 
		    (ftnlen)6, (ftnlen)2);
#line 251 "cunmql.f"
	    nb = min(i__1,i__2);
#line 253 "cunmql.f"
	    lwkopt = nw * nb + 4160;
#line 254 "cunmql.f"
	}
#line 255 "cunmql.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 256 "cunmql.f"
    }

#line 258 "cunmql.f"
    if (*info != 0) {
#line 259 "cunmql.f"
	i__1 = -(*info);
#line 259 "cunmql.f"
	xerbla_("CUNMQL", &i__1, (ftnlen)6);
#line 260 "cunmql.f"
	return 0;
#line 261 "cunmql.f"
    } else if (lquery) {
#line 262 "cunmql.f"
	return 0;
#line 263 "cunmql.f"
    }

/*     Quick return if possible */

#line 267 "cunmql.f"
    if (*m == 0 || *n == 0) {
#line 268 "cunmql.f"
	return 0;
#line 269 "cunmql.f"
    }

/*     Determine the block size */

#line 273 "cunmql.f"
    nbmin = 2;
#line 274 "cunmql.f"
    ldwork = nw;
#line 275 "cunmql.f"
    if (nb > 1 && nb < *k) {
#line 276 "cunmql.f"
	if (*lwork < nw * nb + 4160) {
#line 277 "cunmql.f"
	    nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
#line 278 "cunmql.f"
	    i__3[0] = 1, a__1[0] = side;
#line 278 "cunmql.f"
	    i__3[1] = 1, a__1[1] = trans;
#line 278 "cunmql.f"
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 278 "cunmql.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CUNMQL", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 278 "cunmql.f"
	    nbmin = max(i__1,i__2);
#line 280 "cunmql.f"
	}
#line 281 "cunmql.f"
    }

#line 283 "cunmql.f"
    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

#line 287 "cunmql.f"
	cunm2l_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
#line 289 "cunmql.f"
    } else {

/*        Use blocked code */

#line 293 "cunmql.f"
	iwt = nw * nb + 1;
#line 294 "cunmql.f"
	if (left && notran || ! left && ! notran) {
#line 296 "cunmql.f"
	    i1 = 1;
#line 297 "cunmql.f"
	    i2 = *k;
#line 298 "cunmql.f"
	    i3 = nb;
#line 299 "cunmql.f"
	} else {
#line 300 "cunmql.f"
	    i1 = (*k - 1) / nb * nb + 1;
#line 301 "cunmql.f"
	    i2 = 1;
#line 302 "cunmql.f"
	    i3 = -nb;
#line 303 "cunmql.f"
	}

#line 305 "cunmql.f"
	if (left) {
#line 306 "cunmql.f"
	    ni = *n;
#line 307 "cunmql.f"
	} else {
#line 308 "cunmql.f"
	    mi = *m;
#line 309 "cunmql.f"
	}

#line 311 "cunmql.f"
	i__1 = i2;
#line 311 "cunmql.f"
	i__2 = i3;
#line 311 "cunmql.f"
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 312 "cunmql.f"
	    i__4 = nb, i__5 = *k - i__ + 1;
#line 312 "cunmql.f"
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

#line 317 "cunmql.f"
	    i__4 = nq - *k + i__ + ib - 1;
#line 317 "cunmql.f"
	    clarft_("Backward", "Columnwise", &i__4, &ib, &a[i__ * a_dim1 + 1]
		    , lda, &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)
		    10);
#line 319 "cunmql.f"
	    if (left) {

/*              H or H**H is applied to C(1:m-k+i+ib-1,1:n) */

#line 323 "cunmql.f"
		mi = *m - *k + i__ + ib - 1;
#line 324 "cunmql.f"
	    } else {

/*              H or H**H is applied to C(1:m,1:n-k+i+ib-1) */

#line 328 "cunmql.f"
		ni = *n - *k + i__ + ib - 1;
#line 329 "cunmql.f"
	    }

/*           Apply H or H**H */

#line 333 "cunmql.f"
	    clarfb_(side, trans, "Backward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ * a_dim1 + 1], lda, &work[iwt], &c__65, &c__[c_offset]
		    , ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8,
		     (ftnlen)10);
#line 336 "cunmql.f"
/* L10: */
#line 336 "cunmql.f"
	}
#line 337 "cunmql.f"
    }
#line 338 "cunmql.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 339 "cunmql.f"
    return 0;

/*     End of CUNMQL */

} /* cunmql_ */

