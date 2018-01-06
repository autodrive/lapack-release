#line 1 "sormbr.f"
/* sormbr.f -- translated by f2c (version 20100827).
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

#line 1 "sormbr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b SORMBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORMBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, */
/*                          LDC, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, VECT */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > If VECT = 'Q', SORMBR overwrites the general real M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > If VECT = 'P', SORMBR overwrites the general real M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      P * C          C * P */
/* > TRANS = 'T':      P**T * C       C * P**T */
/* > */
/* > Here Q and P**T are the orthogonal matrices determined by SGEBRD when */
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
/* >          matrix reduced by SGEBRD. */
/* >          If VECT = 'P', the number of rows in the original */
/* >          matrix reduced by SGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension */
/* >                                (LDA,min(nq,K)) if VECT = 'Q' */
/* >                                (LDA,nq)        if VECT = 'P' */
/* >          The vectors which define the elementary reflectors H(i) and */
/* >          G(i), whose products determine the matrices Q and P, as */
/* >          returned by SGEBRD. */
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
/* >          TAU is REAL array, dimension (min(nq,K)) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i) which determines Q or P, as returned */
/* >          by SGEBRD in the array argument TAUQ or TAUP. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
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

/*  ===================================================================== */
/* Subroutine */ int sormbr_(char *vect, char *side, char *trans, integer *m, 
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
    static logical notran, applyq;
    static char transt[1];
    extern /* Subroutine */ int sormlq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 235 "sormbr.f"
    /* Parameter adjustments */
#line 235 "sormbr.f"
    a_dim1 = *lda;
#line 235 "sormbr.f"
    a_offset = 1 + a_dim1;
#line 235 "sormbr.f"
    a -= a_offset;
#line 235 "sormbr.f"
    --tau;
#line 235 "sormbr.f"
    c_dim1 = *ldc;
#line 235 "sormbr.f"
    c_offset = 1 + c_dim1;
#line 235 "sormbr.f"
    c__ -= c_offset;
#line 235 "sormbr.f"
    --work;
#line 235 "sormbr.f"

#line 235 "sormbr.f"
    /* Function Body */
#line 235 "sormbr.f"
    *info = 0;
#line 236 "sormbr.f"
    applyq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 237 "sormbr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 238 "sormbr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 239 "sormbr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q or P and NW is the minimum dimension of WORK */

#line 243 "sormbr.f"
    if (left) {
#line 244 "sormbr.f"
	nq = *m;
#line 245 "sormbr.f"
	nw = *n;
#line 246 "sormbr.f"
    } else {
#line 247 "sormbr.f"
	nq = *n;
#line 248 "sormbr.f"
	nw = *m;
#line 249 "sormbr.f"
    }
#line 250 "sormbr.f"
    if (! applyq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 251 "sormbr.f"
	*info = -1;
#line 252 "sormbr.f"
    } else if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 253 "sormbr.f"
	*info = -2;
#line 254 "sormbr.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 255 "sormbr.f"
	*info = -3;
#line 256 "sormbr.f"
    } else if (*m < 0) {
#line 257 "sormbr.f"
	*info = -4;
#line 258 "sormbr.f"
    } else if (*n < 0) {
#line 259 "sormbr.f"
	*info = -5;
#line 260 "sormbr.f"
    } else if (*k < 0) {
#line 261 "sormbr.f"
	*info = -6;
#line 262 "sormbr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 262 "sormbr.f"
	i__1 = 1, i__2 = min(nq,*k);
#line 262 "sormbr.f"
	if (applyq && *lda < max(1,nq) || ! applyq && *lda < max(i__1,i__2)) {
#line 265 "sormbr.f"
	    *info = -8;
#line 266 "sormbr.f"
	} else if (*ldc < max(1,*m)) {
#line 267 "sormbr.f"
	    *info = -11;
#line 268 "sormbr.f"
	} else if (*lwork < max(1,nw) && ! lquery) {
#line 269 "sormbr.f"
	    *info = -13;
#line 270 "sormbr.f"
	}
#line 270 "sormbr.f"
    }

#line 272 "sormbr.f"
    if (*info == 0) {
#line 273 "sormbr.f"
	if (applyq) {
#line 274 "sormbr.f"
	    if (left) {
/* Writing concatenation */
#line 275 "sormbr.f"
		i__3[0] = 1, a__1[0] = side;
#line 275 "sormbr.f"
		i__3[1] = 1, a__1[1] = trans;
#line 275 "sormbr.f"
		s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 275 "sormbr.f"
		i__1 = *m - 1;
#line 275 "sormbr.f"
		i__2 = *m - 1;
#line 275 "sormbr.f"
		nb = ilaenv_(&c__1, "SORMQR", ch__1, &i__1, n, &i__2, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 277 "sormbr.f"
	    } else {
/* Writing concatenation */
#line 278 "sormbr.f"
		i__3[0] = 1, a__1[0] = side;
#line 278 "sormbr.f"
		i__3[1] = 1, a__1[1] = trans;
#line 278 "sormbr.f"
		s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 278 "sormbr.f"
		i__1 = *n - 1;
#line 278 "sormbr.f"
		i__2 = *n - 1;
#line 278 "sormbr.f"
		nb = ilaenv_(&c__1, "SORMQR", ch__1, m, &i__1, &i__2, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 280 "sormbr.f"
	    }
#line 281 "sormbr.f"
	} else {
#line 282 "sormbr.f"
	    if (left) {
/* Writing concatenation */
#line 283 "sormbr.f"
		i__3[0] = 1, a__1[0] = side;
#line 283 "sormbr.f"
		i__3[1] = 1, a__1[1] = trans;
#line 283 "sormbr.f"
		s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 283 "sormbr.f"
		i__1 = *m - 1;
#line 283 "sormbr.f"
		i__2 = *m - 1;
#line 283 "sormbr.f"
		nb = ilaenv_(&c__1, "SORMLQ", ch__1, &i__1, n, &i__2, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 285 "sormbr.f"
	    } else {
/* Writing concatenation */
#line 286 "sormbr.f"
		i__3[0] = 1, a__1[0] = side;
#line 286 "sormbr.f"
		i__3[1] = 1, a__1[1] = trans;
#line 286 "sormbr.f"
		s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 286 "sormbr.f"
		i__1 = *n - 1;
#line 286 "sormbr.f"
		i__2 = *n - 1;
#line 286 "sormbr.f"
		nb = ilaenv_(&c__1, "SORMLQ", ch__1, m, &i__1, &i__2, &c_n1, (
			ftnlen)6, (ftnlen)2);
#line 288 "sormbr.f"
	    }
#line 289 "sormbr.f"
	}
#line 290 "sormbr.f"
	lwkopt = max(1,nw) * nb;
#line 291 "sormbr.f"
	work[1] = (doublereal) lwkopt;
#line 292 "sormbr.f"
    }

#line 294 "sormbr.f"
    if (*info != 0) {
#line 295 "sormbr.f"
	i__1 = -(*info);
#line 295 "sormbr.f"
	xerbla_("SORMBR", &i__1, (ftnlen)6);
#line 296 "sormbr.f"
	return 0;
#line 297 "sormbr.f"
    } else if (lquery) {
#line 298 "sormbr.f"
	return 0;
#line 299 "sormbr.f"
    }

/*     Quick return if possible */

#line 303 "sormbr.f"
    work[1] = 1.;
#line 304 "sormbr.f"
    if (*m == 0 || *n == 0) {
#line 304 "sormbr.f"
	return 0;
#line 304 "sormbr.f"
    }

#line 307 "sormbr.f"
    if (applyq) {

/*        Apply Q */

#line 311 "sormbr.f"
	if (nq >= *k) {

/*           Q was determined by a call to SGEBRD with nq >= k */

#line 315 "sormbr.f"
	    sormqr_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
		    ftnlen)1);
#line 317 "sormbr.f"
	} else if (nq > 1) {

/*           Q was determined by a call to SGEBRD with nq < k */

#line 321 "sormbr.f"
	    if (left) {
#line 322 "sormbr.f"
		mi = *m - 1;
#line 323 "sormbr.f"
		ni = *n;
#line 324 "sormbr.f"
		i1 = 2;
#line 325 "sormbr.f"
		i2 = 1;
#line 326 "sormbr.f"
	    } else {
#line 327 "sormbr.f"
		mi = *m;
#line 328 "sormbr.f"
		ni = *n - 1;
#line 329 "sormbr.f"
		i1 = 1;
#line 330 "sormbr.f"
		i2 = 2;
#line 331 "sormbr.f"
	    }
#line 332 "sormbr.f"
	    i__1 = nq - 1;
#line 332 "sormbr.f"
	    sormqr_(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1]
		    , &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (
		    ftnlen)1, (ftnlen)1);
#line 334 "sormbr.f"
	}
#line 335 "sormbr.f"
    } else {

/*        Apply P */

#line 339 "sormbr.f"
	if (notran) {
#line 340 "sormbr.f"
	    *(unsigned char *)transt = 'T';
#line 341 "sormbr.f"
	} else {
#line 342 "sormbr.f"
	    *(unsigned char *)transt = 'N';
#line 343 "sormbr.f"
	}
#line 344 "sormbr.f"
	if (nq > *k) {

/*           P was determined by a call to SGEBRD with nq > k */

#line 348 "sormbr.f"
	    sormlq_(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
		    ftnlen)1);
#line 350 "sormbr.f"
	} else if (nq > 1) {

/*           P was determined by a call to SGEBRD with nq <= k */

#line 354 "sormbr.f"
	    if (left) {
#line 355 "sormbr.f"
		mi = *m - 1;
#line 356 "sormbr.f"
		ni = *n;
#line 357 "sormbr.f"
		i1 = 2;
#line 358 "sormbr.f"
		i2 = 1;
#line 359 "sormbr.f"
	    } else {
#line 360 "sormbr.f"
		mi = *m;
#line 361 "sormbr.f"
		ni = *n - 1;
#line 362 "sormbr.f"
		i1 = 1;
#line 363 "sormbr.f"
		i2 = 2;
#line 364 "sormbr.f"
	    }
#line 365 "sormbr.f"
	    i__1 = nq - 1;
#line 365 "sormbr.f"
	    sormlq_(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda,
		     &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &
		    iinfo, (ftnlen)1, (ftnlen)1);
#line 367 "sormbr.f"
	}
#line 368 "sormbr.f"
    }
#line 369 "sormbr.f"
    work[1] = (doublereal) lwkopt;
#line 370 "sormbr.f"
    return 0;

/*     End of SORMBR */

} /* sormbr_ */

