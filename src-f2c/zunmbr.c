#line 1 "zunmbr.f"
/* zunmbr.f -- translated by f2c (version 20100827).
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

#line 1 "zunmbr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b ZUNMBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNMBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, */
/*                          LDC, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, VECT */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > If VECT = 'Q', ZUNMBR overwrites the general complex M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > If VECT = 'P', ZUNMBR overwrites the general complex M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      P * C          C * P */
/* > TRANS = 'C':      P**H * C       C * P**H */
/* > */
/* > Here Q and P**H are the unitary matrices determined by ZGEBRD when */
/* > reducing a complex matrix A to bidiagonal form: A = Q * B * P**H. Q */
/* > and P**H are defined as products of elementary reflectors H(i) and */
/* > G(i) respectively. */
/* > */
/* > Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the */
/* > order of the unitary matrix Q or P**H that is applied. */
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
/* >          = 'Q': apply Q or Q**H; */
/* >          = 'P': apply P or P**H. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q, Q**H, P or P**H from the Left; */
/* >          = 'R': apply Q, Q**H, P or P**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q or P; */
/* >          = 'C':  Conjugate transpose, apply Q**H or P**H. */
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
/* >          matrix reduced by ZGEBRD. */
/* >          If VECT = 'P', the number of rows in the original */
/* >          matrix reduced by ZGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension */
/* >                                (LDA,min(nq,K)) if VECT = 'Q' */
/* >                                (LDA,nq)        if VECT = 'P' */
/* >          The vectors which define the elementary reflectors H(i) and */
/* >          G(i), whose products determine the matrices Q and P, as */
/* >          returned by ZGEBRD. */
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
/* >          TAU is COMPLEX*16 array, dimension (min(nq,K)) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i) which determines Q or P, as returned */
/* >          by ZGEBRD in the array argument TAUQ or TAUP. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q */
/* >          or P*C or P**H*C or C*P or C*P**H. */
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
/* >          if SIDE = 'R', LWORK >= max(1,M); */
/* >          if N = 0 or M = 0, LWORK >= 1. */
/* >          For optimum performance LWORK >= max(1,N*NB) if SIDE = 'L', */
/* >          and LWORK >= max(1,M*NB) if SIDE = 'R', where NB is the */
/* >          optimal blocksize. (NB = 0 if M = 0 or N = 0.) */
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

/*  ===================================================================== */
/* Subroutine */ int zunmbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	lwork, integer *info, ftnlen vect_len, ftnlen side_len, ftnlen 
	trans_len)
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
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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

#line 234 "zunmbr.f"
    /* Parameter adjustments */
#line 234 "zunmbr.f"
    a_dim1 = *lda;
#line 234 "zunmbr.f"
    a_offset = 1 + a_dim1;
#line 234 "zunmbr.f"
    a -= a_offset;
#line 234 "zunmbr.f"
    --tau;
#line 234 "zunmbr.f"
    c_dim1 = *ldc;
#line 234 "zunmbr.f"
    c_offset = 1 + c_dim1;
#line 234 "zunmbr.f"
    c__ -= c_offset;
#line 234 "zunmbr.f"
    --work;
#line 234 "zunmbr.f"

#line 234 "zunmbr.f"
    /* Function Body */
#line 234 "zunmbr.f"
    *info = 0;
#line 235 "zunmbr.f"
    applyq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 236 "zunmbr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 237 "zunmbr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 238 "zunmbr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q or P and NW is the minimum dimension of WORK */

#line 242 "zunmbr.f"
    if (left) {
#line 243 "zunmbr.f"
	nq = *m;
#line 244 "zunmbr.f"
	nw = *n;
#line 245 "zunmbr.f"
    } else {
#line 246 "zunmbr.f"
	nq = *n;
#line 247 "zunmbr.f"
	nw = *m;
#line 248 "zunmbr.f"
    }
#line 249 "zunmbr.f"
    if (*m == 0 || *n == 0) {
#line 250 "zunmbr.f"
	nw = 0;
#line 251 "zunmbr.f"
    }
#line 252 "zunmbr.f"
    if (! applyq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 253 "zunmbr.f"
	*info = -1;
#line 254 "zunmbr.f"
    } else if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 255 "zunmbr.f"
	*info = -2;
#line 256 "zunmbr.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 257 "zunmbr.f"
	*info = -3;
#line 258 "zunmbr.f"
    } else if (*m < 0) {
#line 259 "zunmbr.f"
	*info = -4;
#line 260 "zunmbr.f"
    } else if (*n < 0) {
#line 261 "zunmbr.f"
	*info = -5;
#line 262 "zunmbr.f"
    } else if (*k < 0) {
#line 263 "zunmbr.f"
	*info = -6;
#line 264 "zunmbr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 264 "zunmbr.f"
	i__1 = 1, i__2 = min(nq,*k);
#line 264 "zunmbr.f"
	if (applyq && *lda < max(1,nq) || ! applyq && *lda < max(i__1,i__2)) {
#line 267 "zunmbr.f"
	    *info = -8;
#line 268 "zunmbr.f"
	} else if (*ldc < max(1,*m)) {
#line 269 "zunmbr.f"
	    *info = -11;
#line 270 "zunmbr.f"
	} else if (*lwork < max(1,nw) && ! lquery) {
#line 271 "zunmbr.f"
	    *info = -13;
#line 272 "zunmbr.f"
	}
#line 272 "zunmbr.f"
    }

#line 274 "zunmbr.f"
    if (*info == 0) {
#line 275 "zunmbr.f"
	if (nw > 0) {
#line 276 "zunmbr.f"
	    if (applyq) {
#line 277 "zunmbr.f"
		if (left) {
/* Writing concatenation */
#line 278 "zunmbr.f"
		    i__3[0] = 1, a__1[0] = side;
#line 278 "zunmbr.f"
		    i__3[1] = 1, a__1[1] = trans;
#line 278 "zunmbr.f"
		    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 278 "zunmbr.f"
		    i__1 = *m - 1;
#line 278 "zunmbr.f"
		    i__2 = *m - 1;
#line 278 "zunmbr.f"
		    nb = ilaenv_(&c__1, "ZUNMQR", ch__1, &i__1, n, &i__2, &
			    c_n1, (ftnlen)6, (ftnlen)2);
#line 280 "zunmbr.f"
		} else {
/* Writing concatenation */
#line 281 "zunmbr.f"
		    i__3[0] = 1, a__1[0] = side;
#line 281 "zunmbr.f"
		    i__3[1] = 1, a__1[1] = trans;
#line 281 "zunmbr.f"
		    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 281 "zunmbr.f"
		    i__1 = *n - 1;
#line 281 "zunmbr.f"
		    i__2 = *n - 1;
#line 281 "zunmbr.f"
		    nb = ilaenv_(&c__1, "ZUNMQR", ch__1, m, &i__1, &i__2, &
			    c_n1, (ftnlen)6, (ftnlen)2);
#line 283 "zunmbr.f"
		}
#line 284 "zunmbr.f"
	    } else {
#line 285 "zunmbr.f"
		if (left) {
/* Writing concatenation */
#line 286 "zunmbr.f"
		    i__3[0] = 1, a__1[0] = side;
#line 286 "zunmbr.f"
		    i__3[1] = 1, a__1[1] = trans;
#line 286 "zunmbr.f"
		    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 286 "zunmbr.f"
		    i__1 = *m - 1;
#line 286 "zunmbr.f"
		    i__2 = *m - 1;
#line 286 "zunmbr.f"
		    nb = ilaenv_(&c__1, "ZUNMLQ", ch__1, &i__1, n, &i__2, &
			    c_n1, (ftnlen)6, (ftnlen)2);
#line 288 "zunmbr.f"
		} else {
/* Writing concatenation */
#line 289 "zunmbr.f"
		    i__3[0] = 1, a__1[0] = side;
#line 289 "zunmbr.f"
		    i__3[1] = 1, a__1[1] = trans;
#line 289 "zunmbr.f"
		    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 289 "zunmbr.f"
		    i__1 = *n - 1;
#line 289 "zunmbr.f"
		    i__2 = *n - 1;
#line 289 "zunmbr.f"
		    nb = ilaenv_(&c__1, "ZUNMLQ", ch__1, m, &i__1, &i__2, &
			    c_n1, (ftnlen)6, (ftnlen)2);
#line 291 "zunmbr.f"
		}
#line 292 "zunmbr.f"
	    }
/* Computing MAX */
#line 293 "zunmbr.f"
	    i__1 = 1, i__2 = nw * nb;
#line 293 "zunmbr.f"
	    lwkopt = max(i__1,i__2);
#line 294 "zunmbr.f"
	} else {
#line 295 "zunmbr.f"
	    lwkopt = 1;
#line 296 "zunmbr.f"
	}
#line 297 "zunmbr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 298 "zunmbr.f"
    }

#line 300 "zunmbr.f"
    if (*info != 0) {
#line 301 "zunmbr.f"
	i__1 = -(*info);
#line 301 "zunmbr.f"
	xerbla_("ZUNMBR", &i__1, (ftnlen)6);
#line 302 "zunmbr.f"
	return 0;
#line 303 "zunmbr.f"
    } else if (lquery) {
#line 304 "zunmbr.f"
	return 0;
#line 305 "zunmbr.f"
    }

/*     Quick return if possible */

#line 309 "zunmbr.f"
    if (*m == 0 || *n == 0) {
#line 309 "zunmbr.f"
	return 0;
#line 309 "zunmbr.f"
    }

#line 312 "zunmbr.f"
    if (applyq) {

/*        Apply Q */

#line 316 "zunmbr.f"
	if (nq >= *k) {

/*           Q was determined by a call to ZGEBRD with nq >= k */

#line 320 "zunmbr.f"
	    zunmqr_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
		    ftnlen)1);
#line 322 "zunmbr.f"
	} else if (nq > 1) {

/*           Q was determined by a call to ZGEBRD with nq < k */

#line 326 "zunmbr.f"
	    if (left) {
#line 327 "zunmbr.f"
		mi = *m - 1;
#line 328 "zunmbr.f"
		ni = *n;
#line 329 "zunmbr.f"
		i1 = 2;
#line 330 "zunmbr.f"
		i2 = 1;
#line 331 "zunmbr.f"
	    } else {
#line 332 "zunmbr.f"
		mi = *m;
#line 333 "zunmbr.f"
		ni = *n - 1;
#line 334 "zunmbr.f"
		i1 = 1;
#line 335 "zunmbr.f"
		i2 = 2;
#line 336 "zunmbr.f"
	    }
#line 337 "zunmbr.f"
	    i__1 = nq - 1;
#line 337 "zunmbr.f"
	    zunmqr_(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1]
		    , &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (
		    ftnlen)1, (ftnlen)1);
#line 339 "zunmbr.f"
	}
#line 340 "zunmbr.f"
    } else {

/*        Apply P */

#line 344 "zunmbr.f"
	if (notran) {
#line 345 "zunmbr.f"
	    *(unsigned char *)transt = 'C';
#line 346 "zunmbr.f"
	} else {
#line 347 "zunmbr.f"
	    *(unsigned char *)transt = 'N';
#line 348 "zunmbr.f"
	}
#line 349 "zunmbr.f"
	if (nq > *k) {

/*           P was determined by a call to ZGEBRD with nq > k */

#line 353 "zunmbr.f"
	    zunmlq_(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
		    ftnlen)1);
#line 355 "zunmbr.f"
	} else if (nq > 1) {

/*           P was determined by a call to ZGEBRD with nq <= k */

#line 359 "zunmbr.f"
	    if (left) {
#line 360 "zunmbr.f"
		mi = *m - 1;
#line 361 "zunmbr.f"
		ni = *n;
#line 362 "zunmbr.f"
		i1 = 2;
#line 363 "zunmbr.f"
		i2 = 1;
#line 364 "zunmbr.f"
	    } else {
#line 365 "zunmbr.f"
		mi = *m;
#line 366 "zunmbr.f"
		ni = *n - 1;
#line 367 "zunmbr.f"
		i1 = 1;
#line 368 "zunmbr.f"
		i2 = 2;
#line 369 "zunmbr.f"
	    }
#line 370 "zunmbr.f"
	    i__1 = nq - 1;
#line 370 "zunmbr.f"
	    zunmlq_(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda,
		     &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &
		    iinfo, (ftnlen)1, (ftnlen)1);
#line 372 "zunmbr.f"
	}
#line 373 "zunmbr.f"
    }
#line 374 "zunmbr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 375 "zunmbr.f"
    return 0;

/*     End of ZUNMBR */

} /* zunmbr_ */

