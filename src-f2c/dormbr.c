#line 1 "dormbr.f"
/* dormbr.f -- translated by f2c (version 20100827).
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

#line 1 "dormbr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b DORMBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORMBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, */
/*                          LDC, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, VECT */
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
/* > If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      P * C          C * P */
/* > TRANS = 'T':      P**T * C       C * P**T */
/* > */
/* > Here Q and P**T are the orthogonal matrices determined by DGEBRD when */
/* > reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and */
/* > P**T are defined as products of elementary reflectors H(i) and G(i) */
/* > respectively. */
/* > */
/* > Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the */
/* > order of the orthogonal matrix Q or P**T that is applied. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an NQ-by-K matrix: */
/* > if nq >= k, Q = H(1) H(2) . . . H(k); */
/* > if nq < k, Q = H(1) H(2) . . . H(nq-1). */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-NQ matrix: */
/* > if k < nq, P = G(1) G(2) . . . G(k); */
/* > if k >= nq, P = G(1) G(2) . . . G(nq-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          = 'Q': apply Q or Q**T; */
/* >          = 'P': apply P or P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q, Q**T, P or P**T from the Left; */
/* >          = 'R': apply Q, Q**T, P or P**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q  or P; */
/* >          = 'T':  Transpose, apply Q**T or P**T. */
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
/* >          If VECT = 'Q', the number of columns in the original */
/* >          matrix reduced by DGEBRD. */
/* >          If VECT = 'P', the number of rows in the original */
/* >          matrix reduced by DGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension */
/* >                                (LDA,min(nq,K)) if VECT = 'Q' */
/* >                                (LDA,nq)        if VECT = 'P' */
/* >          The vectors which define the elementary reflectors H(i) and */
/* >          G(i), whose products determine the matrices Q and P, as */
/* >          returned by DGEBRD. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If VECT = 'Q', LDA >= max(1,nq); */
/* >          if VECT = 'P', LDA >= max(1,min(nq,K)). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (min(nq,K)) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i) which determines Q or P, as returned */
/* >          by DGEBRD in the array argument TAUQ or TAUP. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q */
/* >          or P*C or P**T*C or C*P or C*P**T. */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dormbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	integer *info, ftnlen vect_len, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2];
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i1, i2, nb, mi, ni, nq, nw;
    static logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dormlq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical notran;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical applyq;
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

#line 233 "dormbr.f"
    /* Parameter adjustments */
#line 233 "dormbr.f"
    a_dim1 = *lda;
#line 233 "dormbr.f"
    a_offset = 1 + a_dim1;
#line 233 "dormbr.f"
    a -= a_offset;
#line 233 "dormbr.f"
    --tau;
#line 233 "dormbr.f"
    c_dim1 = *ldc;
#line 233 "dormbr.f"
    c_offset = 1 + c_dim1;
#line 233 "dormbr.f"
    c__ -= c_offset;
#line 233 "dormbr.f"
    --work;
#line 233 "dormbr.f"

#line 233 "dormbr.f"
    /* Function Body */
#line 233 "dormbr.f"
    *info = 0;
#line 234 "dormbr.f"
    applyq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 235 "dormbr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 236 "dormbr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 237 "dormbr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q or P and NW is the minimum dimension of WORK */

#line 241 "dormbr.f"
    if (left) {
#line 242 "dormbr.f"
	nq = *m;
#line 243 "dormbr.f"
	nw = *n;
#line 244 "dormbr.f"
    } else {
#line 245 "dormbr.f"
	nq = *n;
#line 246 "dormbr.f"
	nw = *m;
#line 247 "dormbr.f"
    }
#line 248 "dormbr.f"
    if (! applyq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 249 "dormbr.f"
	*info = -1;
#line 250 "dormbr.f"
    } else if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 251 "dormbr.f"
	*info = -2;
#line 252 "dormbr.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 253 "dormbr.f"
	*info = -3;
#line 254 "dormbr.f"
    } else if (*m < 0) {
#line 255 "dormbr.f"
	*info = -4;
#line 256 "dormbr.f"
    } else if (*n < 0) {
#line 257 "dormbr.f"
	*info = -5;
#line 258 "dormbr.f"
    } else if (*k < 0) {
#line 259 "dormbr.f"
	*info = -6;
#line 260 "dormbr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 260 "dormbr.f"
	i__1 = 1, i__2 = min(nq,*k);
#line 260 "dormbr.f"
	if (applyq && *lda < max(1,nq) || ! applyq && *lda < max(i__1,i__2)) {
#line 263 "dormbr.f"
	    *info = -8;
#line 264 "dormbr.f"
	} else if (*ldc < max(1,*m)) {
#line 265 "dormbr.f"
	    *info = -11;
#line 266 "dormbr.f"
	} else if (*lwork < max(1,nw) && ! lquery) {
#line 267 "dormbr.f"
	    *info = -13;
#line 268 "dormbr.f"
	}
#line 268 "dormbr.f"
    }

#line 270 "dormbr.f"
    if (*info == 0) {
#line 271 "dormbr.f"
	if (applyq) {
#line 272 "dormbr.f"
	    if (left) {
/* Writing concatenation */
#line 273 "dormbr.f"
		i__3[0] = 1, a__1[0] = side;
#line 273 "dormbr.f"
		i__3[1] = 1, a__1[1] = trans;
#line 273 "dormbr.f"
		s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 273 "dormbr.f"
		i__1 = *m - 1;
#line 273 "dormbr.f"
		i__2 = *m - 1;
#line 273 "dormbr.f"
		nb = ilaenv_(&c__1, "DORMQR", ch__1, &i__1, n, &i__2, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 275 "dormbr.f"
	    } else {
/* Writing concatenation */
#line 276 "dormbr.f"
		i__3[0] = 1, a__1[0] = side;
#line 276 "dormbr.f"
		i__3[1] = 1, a__1[1] = trans;
#line 276 "dormbr.f"
		s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 276 "dormbr.f"
		i__1 = *n - 1;
#line 276 "dormbr.f"
		i__2 = *n - 1;
#line 276 "dormbr.f"
		nb = ilaenv_(&c__1, "DORMQR", ch__1, m, &i__1, &i__2, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 278 "dormbr.f"
	    }
#line 279 "dormbr.f"
	} else {
#line 280 "dormbr.f"
	    if (left) {
/* Writing concatenation */
#line 281 "dormbr.f"
		i__3[0] = 1, a__1[0] = side;
#line 281 "dormbr.f"
		i__3[1] = 1, a__1[1] = trans;
#line 281 "dormbr.f"
		s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 281 "dormbr.f"
		i__1 = *m - 1;
#line 281 "dormbr.f"
		i__2 = *m - 1;
#line 281 "dormbr.f"
		nb = ilaenv_(&c__1, "DORMLQ", ch__1, &i__1, n, &i__2, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 283 "dormbr.f"
	    } else {
/* Writing concatenation */
#line 284 "dormbr.f"
		i__3[0] = 1, a__1[0] = side;
#line 284 "dormbr.f"
		i__3[1] = 1, a__1[1] = trans;
#line 284 "dormbr.f"
		s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 284 "dormbr.f"
		i__1 = *n - 1;
#line 284 "dormbr.f"
		i__2 = *n - 1;
#line 284 "dormbr.f"
		nb = ilaenv_(&c__1, "DORMLQ", ch__1, m, &i__1, &i__2, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 286 "dormbr.f"
	    }
#line 287 "dormbr.f"
	}
#line 288 "dormbr.f"
	lwkopt = max(1,nw) * nb;
#line 289 "dormbr.f"
	work[1] = (doublereal) lwkopt;
#line 290 "dormbr.f"
    }

#line 292 "dormbr.f"
    if (*info != 0) {
#line 293 "dormbr.f"
	i__1 = -(*info);
#line 293 "dormbr.f"
	xerbla_("DORMBR", &i__1, (ftnlen)6);
#line 294 "dormbr.f"
	return 0;
#line 295 "dormbr.f"
    } else if (lquery) {
#line 296 "dormbr.f"
	return 0;
#line 297 "dormbr.f"
    }

/*     Quick return if possible */

#line 301 "dormbr.f"
    work[1] = 1.;
#line 302 "dormbr.f"
    if (*m == 0 || *n == 0) {
#line 302 "dormbr.f"
	return 0;
#line 302 "dormbr.f"
    }

#line 305 "dormbr.f"
    if (applyq) {

/*        Apply Q */

#line 309 "dormbr.f"
	if (nq >= *k) {

/*           Q was determined by a call to DGEBRD with nq >= k */

#line 313 "dormbr.f"
	    dormqr_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
		    ftnlen)1);
#line 315 "dormbr.f"
	} else if (nq > 1) {

/*           Q was determined by a call to DGEBRD with nq < k */

#line 319 "dormbr.f"
	    if (left) {
#line 320 "dormbr.f"
		mi = *m - 1;
#line 321 "dormbr.f"
		ni = *n;
#line 322 "dormbr.f"
		i1 = 2;
#line 323 "dormbr.f"
		i2 = 1;
#line 324 "dormbr.f"
	    } else {
#line 325 "dormbr.f"
		mi = *m;
#line 326 "dormbr.f"
		ni = *n - 1;
#line 327 "dormbr.f"
		i1 = 1;
#line 328 "dormbr.f"
		i2 = 2;
#line 329 "dormbr.f"
	    }
#line 330 "dormbr.f"
	    i__1 = nq - 1;
#line 330 "dormbr.f"
	    dormqr_(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1]
		    , &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (
		    ftnlen)1, (ftnlen)1);
#line 332 "dormbr.f"
	}
#line 333 "dormbr.f"
    } else {

/*        Apply P */

#line 337 "dormbr.f"
	if (notran) {
#line 338 "dormbr.f"
	    *(unsigned char *)transt = 'T';
#line 339 "dormbr.f"
	} else {
#line 340 "dormbr.f"
	    *(unsigned char *)transt = 'N';
#line 341 "dormbr.f"
	}
#line 342 "dormbr.f"
	if (nq > *k) {

/*           P was determined by a call to DGEBRD with nq > k */

#line 346 "dormbr.f"
	    dormlq_(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
		    ftnlen)1);
#line 348 "dormbr.f"
	} else if (nq > 1) {

/*           P was determined by a call to DGEBRD with nq <= k */

#line 352 "dormbr.f"
	    if (left) {
#line 353 "dormbr.f"
		mi = *m - 1;
#line 354 "dormbr.f"
		ni = *n;
#line 355 "dormbr.f"
		i1 = 2;
#line 356 "dormbr.f"
		i2 = 1;
#line 357 "dormbr.f"
	    } else {
#line 358 "dormbr.f"
		mi = *m;
#line 359 "dormbr.f"
		ni = *n - 1;
#line 360 "dormbr.f"
		i1 = 1;
#line 361 "dormbr.f"
		i2 = 2;
#line 362 "dormbr.f"
	    }
#line 363 "dormbr.f"
	    i__1 = nq - 1;
#line 363 "dormbr.f"
	    dormlq_(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda,
		     &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &
		    iinfo, (ftnlen)1, (ftnlen)1);
#line 365 "dormbr.f"
	}
#line 366 "dormbr.f"
    }
#line 367 "dormbr.f"
    work[1] = (doublereal) lwkopt;
#line 368 "dormbr.f"
    return 0;

/*     End of DORMBR */

} /* dormbr_ */

