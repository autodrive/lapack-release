#line 1 "sorgbr.f"
/* sorgbr.f -- translated by f2c (version 20100827).
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

#line 1 "sorgbr.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief \b SORGBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORGBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          VECT */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORGBR generates one of the real orthogonal matrices Q or P**T */
/* > determined by SGEBRD when reducing a real matrix A to bidiagonal */
/* > form: A = Q * B * P**T.  Q and P**T are defined as products of */
/* > elementary reflectors H(i) or G(i) respectively. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q */
/* > is of order M: */
/* > if m >= k, Q = H(1) H(2) . . . H(k) and SORGBR returns the first n */
/* > columns of Q, where m >= n >= k; */
/* > if m < k, Q = H(1) H(2) . . . H(m-1) and SORGBR returns Q as an */
/* > M-by-M matrix. */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T */
/* > is of order N: */
/* > if k < n, P**T = G(k) . . . G(2) G(1) and SORGBR returns the first m */
/* > rows of P**T, where n >= m >= k; */
/* > if k >= n, P**T = G(n-1) . . . G(2) G(1) and SORGBR returns P**T as */
/* > an N-by-N matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          Specifies whether the matrix Q or the matrix P**T is */
/* >          required, as defined in the transformation applied by SGEBRD: */
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
/* >          matrix reduced by SGEBRD. */
/* >          If VECT = 'P', the number of rows in the original K-by-N */
/* >          matrix reduced by SGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by SGEBRD. */
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
/* >          TAU is REAL array, dimension */
/* >                                (min(M,K)) if VECT = 'Q' */
/* >                                (min(N,K)) if VECT = 'P' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i), which determines Q or P**T, as */
/* >          returned by SGEBRD in its array argument TAUQ or TAUP. */
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

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorgbr_(char *vect, integer *m, integer *n, integer *k, 
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), sorglq_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), sorgqr_(
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

#line 197 "sorgbr.f"
    /* Parameter adjustments */
#line 197 "sorgbr.f"
    a_dim1 = *lda;
#line 197 "sorgbr.f"
    a_offset = 1 + a_dim1;
#line 197 "sorgbr.f"
    a -= a_offset;
#line 197 "sorgbr.f"
    --tau;
#line 197 "sorgbr.f"
    --work;
#line 197 "sorgbr.f"

#line 197 "sorgbr.f"
    /* Function Body */
#line 197 "sorgbr.f"
    *info = 0;
#line 198 "sorgbr.f"
    wantq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 199 "sorgbr.f"
    mn = min(*m,*n);
#line 200 "sorgbr.f"
    lquery = *lwork == -1;
#line 201 "sorgbr.f"
    if (! wantq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 202 "sorgbr.f"
	*info = -1;
#line 203 "sorgbr.f"
    } else if (*m < 0) {
#line 204 "sorgbr.f"
	*info = -2;
#line 205 "sorgbr.f"
    } else if (*n < 0 || wantq && (*n > *m || *n < min(*m,*k)) || ! wantq && (
	    *m > *n || *m < min(*n,*k))) {
#line 208 "sorgbr.f"
	*info = -3;
#line 209 "sorgbr.f"
    } else if (*k < 0) {
#line 210 "sorgbr.f"
	*info = -4;
#line 211 "sorgbr.f"
    } else if (*lda < max(1,*m)) {
#line 212 "sorgbr.f"
	*info = -6;
#line 213 "sorgbr.f"
    } else if (*lwork < max(1,mn) && ! lquery) {
#line 214 "sorgbr.f"
	*info = -9;
#line 215 "sorgbr.f"
    }

#line 217 "sorgbr.f"
    if (*info == 0) {
#line 218 "sorgbr.f"
	work[1] = 1.;
#line 219 "sorgbr.f"
	if (wantq) {
#line 220 "sorgbr.f"
	    if (*m >= *k) {
#line 221 "sorgbr.f"
		sorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 222 "sorgbr.f"
	    } else {
#line 223 "sorgbr.f"
		if (*m > 1) {
#line 224 "sorgbr.f"
		    i__1 = *m - 1;
#line 224 "sorgbr.f"
		    i__2 = *m - 1;
#line 224 "sorgbr.f"
		    i__3 = *m - 1;
#line 224 "sorgbr.f"
		    sorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 226 "sorgbr.f"
		}
#line 227 "sorgbr.f"
	    }
#line 228 "sorgbr.f"
	} else {
#line 229 "sorgbr.f"
	    if (*k < *n) {
#line 230 "sorgbr.f"
		sorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 231 "sorgbr.f"
	    } else {
#line 232 "sorgbr.f"
		if (*n > 1) {
#line 233 "sorgbr.f"
		    i__1 = *n - 1;
#line 233 "sorgbr.f"
		    i__2 = *n - 1;
#line 233 "sorgbr.f"
		    i__3 = *n - 1;
#line 233 "sorgbr.f"
		    sorglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 235 "sorgbr.f"
		}
#line 236 "sorgbr.f"
	    }
#line 237 "sorgbr.f"
	}
#line 238 "sorgbr.f"
	lwkopt = (integer) work[1];
#line 239 "sorgbr.f"
	lwkopt = max(lwkopt,mn);
#line 240 "sorgbr.f"
    }

#line 242 "sorgbr.f"
    if (*info != 0) {
#line 243 "sorgbr.f"
	i__1 = -(*info);
#line 243 "sorgbr.f"
	xerbla_("SORGBR", &i__1, (ftnlen)6);
#line 244 "sorgbr.f"
	return 0;
#line 245 "sorgbr.f"
    } else if (lquery) {
#line 246 "sorgbr.f"
	work[1] = (doublereal) lwkopt;
#line 247 "sorgbr.f"
	return 0;
#line 248 "sorgbr.f"
    }

/*     Quick return if possible */

#line 252 "sorgbr.f"
    if (*m == 0 || *n == 0) {
#line 253 "sorgbr.f"
	work[1] = 1.;
#line 254 "sorgbr.f"
	return 0;
#line 255 "sorgbr.f"
    }

#line 257 "sorgbr.f"
    if (wantq) {

/*        Form Q, determined by a call to SGEBRD to reduce an m-by-k */
/*        matrix */

#line 262 "sorgbr.f"
	if (*m >= *k) {

/*           If m >= k, assume m >= n >= k */

#line 266 "sorgbr.f"
	    sorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 268 "sorgbr.f"
	} else {

/*           If m < k, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           column to the right, and set the first row and column of Q */
/*           to those of the unit matrix */

#line 276 "sorgbr.f"
	    for (j = *m; j >= 2; --j) {
#line 277 "sorgbr.f"
		a[j * a_dim1 + 1] = 0.;
#line 278 "sorgbr.f"
		i__1 = *m;
#line 278 "sorgbr.f"
		for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 279 "sorgbr.f"
		    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
#line 280 "sorgbr.f"
/* L10: */
#line 280 "sorgbr.f"
		}
#line 281 "sorgbr.f"
/* L20: */
#line 281 "sorgbr.f"
	    }
#line 282 "sorgbr.f"
	    a[a_dim1 + 1] = 1.;
#line 283 "sorgbr.f"
	    i__1 = *m;
#line 283 "sorgbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 284 "sorgbr.f"
		a[i__ + a_dim1] = 0.;
#line 285 "sorgbr.f"
/* L30: */
#line 285 "sorgbr.f"
	    }
#line 286 "sorgbr.f"
	    if (*m > 1) {

/*              Form Q(2:m,2:m) */

#line 290 "sorgbr.f"
		i__1 = *m - 1;
#line 290 "sorgbr.f"
		i__2 = *m - 1;
#line 290 "sorgbr.f"
		i__3 = *m - 1;
#line 290 "sorgbr.f"
		sorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 292 "sorgbr.f"
	    }
#line 293 "sorgbr.f"
	}
#line 294 "sorgbr.f"
    } else {

/*        Form P**T, determined by a call to SGEBRD to reduce a k-by-n */
/*        matrix */

#line 299 "sorgbr.f"
	if (*k < *n) {

/*           If k < n, assume k <= m <= n */

#line 303 "sorgbr.f"
	    sorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 305 "sorgbr.f"
	} else {

/*           If k >= n, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           row downward, and set the first row and column of P**T to */
/*           those of the unit matrix */

#line 313 "sorgbr.f"
	    a[a_dim1 + 1] = 1.;
#line 314 "sorgbr.f"
	    i__1 = *n;
#line 314 "sorgbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 315 "sorgbr.f"
		a[i__ + a_dim1] = 0.;
#line 316 "sorgbr.f"
/* L40: */
#line 316 "sorgbr.f"
	    }
#line 317 "sorgbr.f"
	    i__1 = *n;
#line 317 "sorgbr.f"
	    for (j = 2; j <= i__1; ++j) {
#line 318 "sorgbr.f"
		for (i__ = j - 1; i__ >= 2; --i__) {
#line 319 "sorgbr.f"
		    a[i__ + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 320 "sorgbr.f"
/* L50: */
#line 320 "sorgbr.f"
		}
#line 321 "sorgbr.f"
		a[j * a_dim1 + 1] = 0.;
#line 322 "sorgbr.f"
/* L60: */
#line 322 "sorgbr.f"
	    }
#line 323 "sorgbr.f"
	    if (*n > 1) {

/*              Form P**T(2:n,2:n) */

#line 327 "sorgbr.f"
		i__1 = *n - 1;
#line 327 "sorgbr.f"
		i__2 = *n - 1;
#line 327 "sorgbr.f"
		i__3 = *n - 1;
#line 327 "sorgbr.f"
		sorglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 329 "sorgbr.f"
	    }
#line 330 "sorgbr.f"
	}
#line 331 "sorgbr.f"
    }
#line 332 "sorgbr.f"
    work[1] = (doublereal) lwkopt;
#line 333 "sorgbr.f"
    return 0;

/*     End of SORGBR */

} /* sorgbr_ */

