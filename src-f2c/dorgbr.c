#line 1 "dorgbr.f"
/* dorgbr.f -- translated by f2c (version 20100827).
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

#line 1 "dorgbr.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief \b DORGBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORGBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          VECT */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORGBR generates one of the real orthogonal matrices Q or P**T */
/* > determined by DGEBRD when reducing a real matrix A to bidiagonal */
/* > form: A = Q * B * P**T.  Q and P**T are defined as products of */
/* > elementary reflectors H(i) or G(i) respectively. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q */
/* > is of order M: */
/* > if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n */
/* > columns of Q, where m >= n >= k; */
/* > if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an */
/* > M-by-M matrix. */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T */
/* > is of order N: */
/* > if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m */
/* > rows of P**T, where n >= m >= k; */
/* > if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as */
/* > an N-by-N matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          Specifies whether the matrix Q or the matrix P**T is */
/* >          required, as defined in the transformation applied by DGEBRD: */
/* >          = 'Q':  generate Q; */
/* >          = 'P':  generate P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix Q or P**T to be returned. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix Q or P**T to be returned. */
/* >          N >= 0. */
/* >          If VECT = 'Q', M >= N >= min(M,K); */
/* >          if VECT = 'P', N >= M >= min(N,K). */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          If VECT = 'Q', the number of columns in the original M-by-K */
/* >          matrix reduced by DGEBRD. */
/* >          If VECT = 'P', the number of rows in the original K-by-N */
/* >          matrix reduced by DGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by DGEBRD. */
/* >          On exit, the M-by-N matrix Q or P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension */
/* >                                (min(M,K)) if VECT = 'Q' */
/* >                                (min(N,K)) if VECT = 'P' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i), which determines Q or P**T, as */
/* >          returned by DGEBRD in its array argument TAUQ or TAUP. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,min(M,N)). */
/* >          For optimum performance LWORK >= min(M,N)*NB, where NB */
/* >          is the optimal blocksize. */
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

/* > \date April 2012 */

/* > \ingroup doubleGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorgbr_(char *vect, integer *m, integer *n, integer *k, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info, ftnlen vect_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, mn;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical wantq;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dorglq_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

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

#line 197 "dorgbr.f"
    /* Parameter adjustments */
#line 197 "dorgbr.f"
    a_dim1 = *lda;
#line 197 "dorgbr.f"
    a_offset = 1 + a_dim1;
#line 197 "dorgbr.f"
    a -= a_offset;
#line 197 "dorgbr.f"
    --tau;
#line 197 "dorgbr.f"
    --work;
#line 197 "dorgbr.f"

#line 197 "dorgbr.f"
    /* Function Body */
#line 197 "dorgbr.f"
    *info = 0;
#line 198 "dorgbr.f"
    wantq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 199 "dorgbr.f"
    mn = min(*m,*n);
#line 200 "dorgbr.f"
    lquery = *lwork == -1;
#line 201 "dorgbr.f"
    if (! wantq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 202 "dorgbr.f"
	*info = -1;
#line 203 "dorgbr.f"
    } else if (*m < 0) {
#line 204 "dorgbr.f"
	*info = -2;
#line 205 "dorgbr.f"
    } else if (*n < 0 || wantq && (*n > *m || *n < min(*m,*k)) || ! wantq && (
	    *m > *n || *m < min(*n,*k))) {
#line 208 "dorgbr.f"
	*info = -3;
#line 209 "dorgbr.f"
    } else if (*k < 0) {
#line 210 "dorgbr.f"
	*info = -4;
#line 211 "dorgbr.f"
    } else if (*lda < max(1,*m)) {
#line 212 "dorgbr.f"
	*info = -6;
#line 213 "dorgbr.f"
    } else if (*lwork < max(1,mn) && ! lquery) {
#line 214 "dorgbr.f"
	*info = -9;
#line 215 "dorgbr.f"
    }

#line 217 "dorgbr.f"
    if (*info == 0) {
#line 218 "dorgbr.f"
	work[1] = 1.;
#line 219 "dorgbr.f"
	if (wantq) {
#line 220 "dorgbr.f"
	    if (*m >= *k) {
#line 221 "dorgbr.f"
		dorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 222 "dorgbr.f"
	    } else {
#line 223 "dorgbr.f"
		if (*m > 1) {
#line 224 "dorgbr.f"
		    i__1 = *m - 1;
#line 224 "dorgbr.f"
		    i__2 = *m - 1;
#line 224 "dorgbr.f"
		    i__3 = *m - 1;
#line 224 "dorgbr.f"
		    dorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 226 "dorgbr.f"
		}
#line 227 "dorgbr.f"
	    }
#line 228 "dorgbr.f"
	} else {
#line 229 "dorgbr.f"
	    if (*k < *n) {
#line 230 "dorgbr.f"
		dorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 231 "dorgbr.f"
	    } else {
#line 232 "dorgbr.f"
		if (*n > 1) {
#line 233 "dorgbr.f"
		    i__1 = *n - 1;
#line 233 "dorgbr.f"
		    i__2 = *n - 1;
#line 233 "dorgbr.f"
		    i__3 = *n - 1;
#line 233 "dorgbr.f"
		    dorglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 235 "dorgbr.f"
		}
#line 236 "dorgbr.f"
	    }
#line 237 "dorgbr.f"
	}
#line 238 "dorgbr.f"
	lwkopt = (integer) work[1];
#line 239 "dorgbr.f"
	lwkopt = max(lwkopt,mn);
#line 240 "dorgbr.f"
    }

#line 242 "dorgbr.f"
    if (*info != 0) {
#line 243 "dorgbr.f"
	i__1 = -(*info);
#line 243 "dorgbr.f"
	xerbla_("DORGBR", &i__1, (ftnlen)6);
#line 244 "dorgbr.f"
	return 0;
#line 245 "dorgbr.f"
    } else if (lquery) {
#line 246 "dorgbr.f"
	work[1] = (doublereal) lwkopt;
#line 247 "dorgbr.f"
	return 0;
#line 248 "dorgbr.f"
    }

/*     Quick return if possible */

#line 252 "dorgbr.f"
    if (*m == 0 || *n == 0) {
#line 253 "dorgbr.f"
	work[1] = 1.;
#line 254 "dorgbr.f"
	return 0;
#line 255 "dorgbr.f"
    }

#line 257 "dorgbr.f"
    if (wantq) {

/*        Form Q, determined by a call to DGEBRD to reduce an m-by-k */
/*        matrix */

#line 262 "dorgbr.f"
	if (*m >= *k) {

/*           If m >= k, assume m >= n >= k */

#line 266 "dorgbr.f"
	    dorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 268 "dorgbr.f"
	} else {

/*           If m < k, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           column to the right, and set the first row and column of Q */
/*           to those of the unit matrix */

#line 276 "dorgbr.f"
	    for (j = *m; j >= 2; --j) {
#line 277 "dorgbr.f"
		a[j * a_dim1 + 1] = 0.;
#line 278 "dorgbr.f"
		i__1 = *m;
#line 278 "dorgbr.f"
		for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 279 "dorgbr.f"
		    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
#line 280 "dorgbr.f"
/* L10: */
#line 280 "dorgbr.f"
		}
#line 281 "dorgbr.f"
/* L20: */
#line 281 "dorgbr.f"
	    }
#line 282 "dorgbr.f"
	    a[a_dim1 + 1] = 1.;
#line 283 "dorgbr.f"
	    i__1 = *m;
#line 283 "dorgbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 284 "dorgbr.f"
		a[i__ + a_dim1] = 0.;
#line 285 "dorgbr.f"
/* L30: */
#line 285 "dorgbr.f"
	    }
#line 286 "dorgbr.f"
	    if (*m > 1) {

/*              Form Q(2:m,2:m) */

#line 290 "dorgbr.f"
		i__1 = *m - 1;
#line 290 "dorgbr.f"
		i__2 = *m - 1;
#line 290 "dorgbr.f"
		i__3 = *m - 1;
#line 290 "dorgbr.f"
		dorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 292 "dorgbr.f"
	    }
#line 293 "dorgbr.f"
	}
#line 294 "dorgbr.f"
    } else {

/*        Form P**T, determined by a call to DGEBRD to reduce a k-by-n */
/*        matrix */

#line 299 "dorgbr.f"
	if (*k < *n) {

/*           If k < n, assume k <= m <= n */

#line 303 "dorgbr.f"
	    dorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 305 "dorgbr.f"
	} else {

/*           If k >= n, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           row downward, and set the first row and column of P**T to */
/*           those of the unit matrix */

#line 313 "dorgbr.f"
	    a[a_dim1 + 1] = 1.;
#line 314 "dorgbr.f"
	    i__1 = *n;
#line 314 "dorgbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 315 "dorgbr.f"
		a[i__ + a_dim1] = 0.;
#line 316 "dorgbr.f"
/* L40: */
#line 316 "dorgbr.f"
	    }
#line 317 "dorgbr.f"
	    i__1 = *n;
#line 317 "dorgbr.f"
	    for (j = 2; j <= i__1; ++j) {
#line 318 "dorgbr.f"
		for (i__ = j - 1; i__ >= 2; --i__) {
#line 319 "dorgbr.f"
		    a[i__ + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 320 "dorgbr.f"
/* L50: */
#line 320 "dorgbr.f"
		}
#line 321 "dorgbr.f"
		a[j * a_dim1 + 1] = 0.;
#line 322 "dorgbr.f"
/* L60: */
#line 322 "dorgbr.f"
	    }
#line 323 "dorgbr.f"
	    if (*n > 1) {

/*              Form P**T(2:n,2:n) */

#line 327 "dorgbr.f"
		i__1 = *n - 1;
#line 327 "dorgbr.f"
		i__2 = *n - 1;
#line 327 "dorgbr.f"
		i__3 = *n - 1;
#line 327 "dorgbr.f"
		dorglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 329 "dorgbr.f"
	    }
#line 330 "dorgbr.f"
	}
#line 331 "dorgbr.f"
    }
#line 332 "dorgbr.f"
    work[1] = (doublereal) lwkopt;
#line 333 "dorgbr.f"
    return 0;

/*     End of DORGBR */

} /* dorgbr_ */

