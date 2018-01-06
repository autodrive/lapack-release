#line 1 "cunmbr.f"
/* cunmbr.f -- translated by f2c (version 20100827).
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

#line 1 "cunmbr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CUNMBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNMBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, */
/*                          LDC, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, VECT */
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
/* > If VECT = 'Q', CUNMBR overwrites the general complex M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'C':      Q**H * C       C * Q**H */
/* > */
/* > If VECT = 'P', CUNMBR overwrites the general complex M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      P * C          C * P */
/* > TRANS = 'C':      P**H * C       C * P**H */
/* > */
/* > Here Q and P**H are the unitary matrices determined by CGEBRD when */
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
/* >          matrix reduced by CGEBRD. */
/* >          If VECT = 'P', the number of rows in the original */
/* >          matrix reduced by CGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension */
/* >                                (LDA,min(nq,K)) if VECT = 'Q' */
/* >                                (LDA,nq)        if VECT = 'P' */
/* >          The vectors which define the elementary reflectors H(i) and */
/* >          G(i), whose products determine the matrices Q and P, as */
/* >          returned by CGEBRD. */
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
/* >          TAU is COMPLEX array, dimension (min(nq,K)) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i) which determines Q or P, as returned */
/* >          by CGEBRD in the array argument TAUQ or TAUP. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC,N) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunmbr_(char *vect, char *side, char *trans, integer *m, 
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
    extern /* Subroutine */ int cunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static logical notran;
    extern /* Subroutine */ int cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static logical applyq;
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

#line 236 "cunmbr.f"
    /* Parameter adjustments */
#line 236 "cunmbr.f"
    a_dim1 = *lda;
#line 236 "cunmbr.f"
    a_offset = 1 + a_dim1;
#line 236 "cunmbr.f"
    a -= a_offset;
#line 236 "cunmbr.f"
    --tau;
#line 236 "cunmbr.f"
    c_dim1 = *ldc;
#line 236 "cunmbr.f"
    c_offset = 1 + c_dim1;
#line 236 "cunmbr.f"
    c__ -= c_offset;
#line 236 "cunmbr.f"
    --work;
#line 236 "cunmbr.f"

#line 236 "cunmbr.f"
    /* Function Body */
#line 236 "cunmbr.f"
    *info = 0;
#line 237 "cunmbr.f"
    applyq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 238 "cunmbr.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 239 "cunmbr.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 240 "cunmbr.f"
    lquery = *lwork == -1;

/*     NQ is the order of Q or P and NW is the minimum dimension of WORK */

#line 244 "cunmbr.f"
    if (left) {
#line 245 "cunmbr.f"
	nq = *m;
#line 246 "cunmbr.f"
	nw = *n;
#line 247 "cunmbr.f"
    } else {
#line 248 "cunmbr.f"
	nq = *n;
#line 249 "cunmbr.f"
	nw = *m;
#line 250 "cunmbr.f"
    }
#line 251 "cunmbr.f"
    if (*m == 0 || *n == 0) {
#line 252 "cunmbr.f"
	nw = 0;
#line 253 "cunmbr.f"
    }
#line 254 "cunmbr.f"
    if (! applyq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 255 "cunmbr.f"
	*info = -1;
#line 256 "cunmbr.f"
    } else if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 257 "cunmbr.f"
	*info = -2;
#line 258 "cunmbr.f"
    } else if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 259 "cunmbr.f"
	*info = -3;
#line 260 "cunmbr.f"
    } else if (*m < 0) {
#line 261 "cunmbr.f"
	*info = -4;
#line 262 "cunmbr.f"
    } else if (*n < 0) {
#line 263 "cunmbr.f"
	*info = -5;
#line 264 "cunmbr.f"
    } else if (*k < 0) {
#line 265 "cunmbr.f"
	*info = -6;
#line 266 "cunmbr.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 266 "cunmbr.f"
	i__1 = 1, i__2 = min(nq,*k);
#line 266 "cunmbr.f"
	if (applyq && *lda < max(1,nq) || ! applyq && *lda < max(i__1,i__2)) {
#line 269 "cunmbr.f"
	    *info = -8;
#line 270 "cunmbr.f"
	} else if (*ldc < max(1,*m)) {
#line 271 "cunmbr.f"
	    *info = -11;
#line 272 "cunmbr.f"
	} else if (*lwork < max(1,nw) && ! lquery) {
#line 273 "cunmbr.f"
	    *info = -13;
#line 274 "cunmbr.f"
	}
#line 274 "cunmbr.f"
    }

#line 276 "cunmbr.f"
    if (*info == 0) {
#line 277 "cunmbr.f"
	if (nw > 0) {
#line 278 "cunmbr.f"
	    if (applyq) {
#line 279 "cunmbr.f"
		if (left) {
/* Writing concatenation */
#line 280 "cunmbr.f"
		    i__3[0] = 1, a__1[0] = side;
#line 280 "cunmbr.f"
		    i__3[1] = 1, a__1[1] = trans;
#line 280 "cunmbr.f"
		    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 280 "cunmbr.f"
		    i__1 = *m - 1;
#line 280 "cunmbr.f"
		    i__2 = *m - 1;
#line 280 "cunmbr.f"
		    nb = ilaenv_(&c__1, "CUNMQR", ch__1, &i__1, n, &i__2, &
			    c_n1, (ftnlen)6, (ftnlen)2);
#line 282 "cunmbr.f"
		} else {
/* Writing concatenation */
#line 283 "cunmbr.f"
		    i__3[0] = 1, a__1[0] = side;
#line 283 "cunmbr.f"
		    i__3[1] = 1, a__1[1] = trans;
#line 283 "cunmbr.f"
		    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 283 "cunmbr.f"
		    i__1 = *n - 1;
#line 283 "cunmbr.f"
		    i__2 = *n - 1;
#line 283 "cunmbr.f"
		    nb = ilaenv_(&c__1, "CUNMQR", ch__1, m, &i__1, &i__2, &
			    c_n1, (ftnlen)6, (ftnlen)2);
#line 285 "cunmbr.f"
		}
#line 286 "cunmbr.f"
	    } else {
#line 287 "cunmbr.f"
		if (left) {
/* Writing concatenation */
#line 288 "cunmbr.f"
		    i__3[0] = 1, a__1[0] = side;
#line 288 "cunmbr.f"
		    i__3[1] = 1, a__1[1] = trans;
#line 288 "cunmbr.f"
		    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 288 "cunmbr.f"
		    i__1 = *m - 1;
#line 288 "cunmbr.f"
		    i__2 = *m - 1;
#line 288 "cunmbr.f"
		    nb = ilaenv_(&c__1, "CUNMLQ", ch__1, &i__1, n, &i__2, &
			    c_n1, (ftnlen)6, (ftnlen)2);
#line 290 "cunmbr.f"
		} else {
/* Writing concatenation */
#line 291 "cunmbr.f"
		    i__3[0] = 1, a__1[0] = side;
#line 291 "cunmbr.f"
		    i__3[1] = 1, a__1[1] = trans;
#line 291 "cunmbr.f"
		    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 291 "cunmbr.f"
		    i__1 = *n - 1;
#line 291 "cunmbr.f"
		    i__2 = *n - 1;
#line 291 "cunmbr.f"
		    nb = ilaenv_(&c__1, "CUNMLQ", ch__1, m, &i__1, &i__2, &
			    c_n1, (ftnlen)6, (ftnlen)2);
#line 293 "cunmbr.f"
		}
#line 294 "cunmbr.f"
	    }
/* Computing MAX */
#line 295 "cunmbr.f"
	    i__1 = 1, i__2 = nw * nb;
#line 295 "cunmbr.f"
	    lwkopt = max(i__1,i__2);
#line 296 "cunmbr.f"
	} else {
#line 297 "cunmbr.f"
	    lwkopt = 1;
#line 298 "cunmbr.f"
	}
#line 299 "cunmbr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 300 "cunmbr.f"
    }

#line 302 "cunmbr.f"
    if (*info != 0) {
#line 303 "cunmbr.f"
	i__1 = -(*info);
#line 303 "cunmbr.f"
	xerbla_("CUNMBR", &i__1, (ftnlen)6);
#line 304 "cunmbr.f"
	return 0;
#line 305 "cunmbr.f"
    } else if (lquery) {
#line 306 "cunmbr.f"
	return 0;
#line 307 "cunmbr.f"
    }

/*     Quick return if possible */

#line 311 "cunmbr.f"
    if (*m == 0 || *n == 0) {
#line 311 "cunmbr.f"
	return 0;
#line 311 "cunmbr.f"
    }

#line 314 "cunmbr.f"
    if (applyq) {

/*        Apply Q */

#line 318 "cunmbr.f"
	if (nq >= *k) {

/*           Q was determined by a call to CGEBRD with nq >= k */

#line 322 "cunmbr.f"
	    cunmqr_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
		    ftnlen)1);
#line 324 "cunmbr.f"
	} else if (nq > 1) {

/*           Q was determined by a call to CGEBRD with nq < k */

#line 328 "cunmbr.f"
	    if (left) {
#line 329 "cunmbr.f"
		mi = *m - 1;
#line 330 "cunmbr.f"
		ni = *n;
#line 331 "cunmbr.f"
		i1 = 2;
#line 332 "cunmbr.f"
		i2 = 1;
#line 333 "cunmbr.f"
	    } else {
#line 334 "cunmbr.f"
		mi = *m;
#line 335 "cunmbr.f"
		ni = *n - 1;
#line 336 "cunmbr.f"
		i1 = 1;
#line 337 "cunmbr.f"
		i2 = 2;
#line 338 "cunmbr.f"
	    }
#line 339 "cunmbr.f"
	    i__1 = nq - 1;
#line 339 "cunmbr.f"
	    cunmqr_(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1]
		    , &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (
		    ftnlen)1, (ftnlen)1);
#line 341 "cunmbr.f"
	}
#line 342 "cunmbr.f"
    } else {

/*        Apply P */

#line 346 "cunmbr.f"
	if (notran) {
#line 347 "cunmbr.f"
	    *(unsigned char *)transt = 'C';
#line 348 "cunmbr.f"
	} else {
#line 349 "cunmbr.f"
	    *(unsigned char *)transt = 'N';
#line 350 "cunmbr.f"
	}
#line 351 "cunmbr.f"
	if (nq > *k) {

/*           P was determined by a call to CGEBRD with nq > k */

#line 355 "cunmbr.f"
	    cunmlq_(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
		    ftnlen)1);
#line 357 "cunmbr.f"
	} else if (nq > 1) {

/*           P was determined by a call to CGEBRD with nq <= k */

#line 361 "cunmbr.f"
	    if (left) {
#line 362 "cunmbr.f"
		mi = *m - 1;
#line 363 "cunmbr.f"
		ni = *n;
#line 364 "cunmbr.f"
		i1 = 2;
#line 365 "cunmbr.f"
		i2 = 1;
#line 366 "cunmbr.f"
	    } else {
#line 367 "cunmbr.f"
		mi = *m;
#line 368 "cunmbr.f"
		ni = *n - 1;
#line 369 "cunmbr.f"
		i1 = 1;
#line 370 "cunmbr.f"
		i2 = 2;
#line 371 "cunmbr.f"
	    }
#line 372 "cunmbr.f"
	    i__1 = nq - 1;
#line 372 "cunmbr.f"
	    cunmlq_(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda,
		     &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &
		    iinfo, (ftnlen)1, (ftnlen)1);
#line 374 "cunmbr.f"
	}
#line 375 "cunmbr.f"
    }
#line 376 "cunmbr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 377 "cunmbr.f"
    return 0;

/*     End of CUNMBR */

} /* cunmbr_ */

