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


/*  -- LAPACK computational routine (version 3.4.1) -- */
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

#line 198 "dorgbr.f"
    /* Parameter adjustments */
#line 198 "dorgbr.f"
    a_dim1 = *lda;
#line 198 "dorgbr.f"
    a_offset = 1 + a_dim1;
#line 198 "dorgbr.f"
    a -= a_offset;
#line 198 "dorgbr.f"
    --tau;
#line 198 "dorgbr.f"
    --work;
#line 198 "dorgbr.f"

#line 198 "dorgbr.f"
    /* Function Body */
#line 198 "dorgbr.f"
    *info = 0;
#line 199 "dorgbr.f"
    wantq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 200 "dorgbr.f"
    mn = min(*m,*n);
#line 201 "dorgbr.f"
    lquery = *lwork == -1;
#line 202 "dorgbr.f"
    if (! wantq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 203 "dorgbr.f"
	*info = -1;
#line 204 "dorgbr.f"
    } else if (*m < 0) {
#line 205 "dorgbr.f"
	*info = -2;
#line 206 "dorgbr.f"
    } else if (*n < 0 || wantq && (*n > *m || *n < min(*m,*k)) || ! wantq && (
	    *m > *n || *m < min(*n,*k))) {
#line 209 "dorgbr.f"
	*info = -3;
#line 210 "dorgbr.f"
    } else if (*k < 0) {
#line 211 "dorgbr.f"
	*info = -4;
#line 212 "dorgbr.f"
    } else if (*lda < max(1,*m)) {
#line 213 "dorgbr.f"
	*info = -6;
#line 214 "dorgbr.f"
    } else if (*lwork < max(1,mn) && ! lquery) {
#line 215 "dorgbr.f"
	*info = -9;
#line 216 "dorgbr.f"
    }

#line 218 "dorgbr.f"
    if (*info == 0) {
#line 219 "dorgbr.f"
	work[1] = 1.;
#line 220 "dorgbr.f"
	if (wantq) {
#line 221 "dorgbr.f"
	    if (*m >= *k) {
#line 222 "dorgbr.f"
		dorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 223 "dorgbr.f"
	    } else {
#line 224 "dorgbr.f"
		if (*m > 1) {
#line 225 "dorgbr.f"
		    i__1 = *m - 1;
#line 225 "dorgbr.f"
		    i__2 = *m - 1;
#line 225 "dorgbr.f"
		    i__3 = *m - 1;
#line 225 "dorgbr.f"
		    dorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 227 "dorgbr.f"
		}
#line 228 "dorgbr.f"
	    }
#line 229 "dorgbr.f"
	} else {
#line 230 "dorgbr.f"
	    if (*k < *n) {
#line 231 "dorgbr.f"
		dorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 232 "dorgbr.f"
	    } else {
#line 233 "dorgbr.f"
		if (*n > 1) {
#line 234 "dorgbr.f"
		    i__1 = *n - 1;
#line 234 "dorgbr.f"
		    i__2 = *n - 1;
#line 234 "dorgbr.f"
		    i__3 = *n - 1;
#line 234 "dorgbr.f"
		    dorglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 236 "dorgbr.f"
		}
#line 237 "dorgbr.f"
	    }
#line 238 "dorgbr.f"
	}
#line 239 "dorgbr.f"
	lwkopt = (integer) work[1];
#line 240 "dorgbr.f"
	lwkopt = max(lwkopt,mn);
#line 241 "dorgbr.f"
    }

#line 243 "dorgbr.f"
    if (*info != 0) {
#line 244 "dorgbr.f"
	i__1 = -(*info);
#line 244 "dorgbr.f"
	xerbla_("DORGBR", &i__1, (ftnlen)6);
#line 245 "dorgbr.f"
	return 0;
#line 246 "dorgbr.f"
    } else if (lquery) {
#line 247 "dorgbr.f"
	work[1] = (doublereal) lwkopt;
#line 248 "dorgbr.f"
	return 0;
#line 249 "dorgbr.f"
    }

/*     Quick return if possible */

#line 253 "dorgbr.f"
    if (*m == 0 || *n == 0) {
#line 254 "dorgbr.f"
	work[1] = 1.;
#line 255 "dorgbr.f"
	return 0;
#line 256 "dorgbr.f"
    }

#line 258 "dorgbr.f"
    if (wantq) {

/*        Form Q, determined by a call to DGEBRD to reduce an m-by-k */
/*        matrix */

#line 263 "dorgbr.f"
	if (*m >= *k) {

/*           If m >= k, assume m >= n >= k */

#line 267 "dorgbr.f"
	    dorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 269 "dorgbr.f"
	} else {

/*           If m < k, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           column to the right, and set the first row and column of Q */
/*           to those of the unit matrix */

#line 277 "dorgbr.f"
	    for (j = *m; j >= 2; --j) {
#line 278 "dorgbr.f"
		a[j * a_dim1 + 1] = 0.;
#line 279 "dorgbr.f"
		i__1 = *m;
#line 279 "dorgbr.f"
		for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 280 "dorgbr.f"
		    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
#line 281 "dorgbr.f"
/* L10: */
#line 281 "dorgbr.f"
		}
#line 282 "dorgbr.f"
/* L20: */
#line 282 "dorgbr.f"
	    }
#line 283 "dorgbr.f"
	    a[a_dim1 + 1] = 1.;
#line 284 "dorgbr.f"
	    i__1 = *m;
#line 284 "dorgbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 285 "dorgbr.f"
		a[i__ + a_dim1] = 0.;
#line 286 "dorgbr.f"
/* L30: */
#line 286 "dorgbr.f"
	    }
#line 287 "dorgbr.f"
	    if (*m > 1) {

/*              Form Q(2:m,2:m) */

#line 291 "dorgbr.f"
		i__1 = *m - 1;
#line 291 "dorgbr.f"
		i__2 = *m - 1;
#line 291 "dorgbr.f"
		i__3 = *m - 1;
#line 291 "dorgbr.f"
		dorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 293 "dorgbr.f"
	    }
#line 294 "dorgbr.f"
	}
#line 295 "dorgbr.f"
    } else {

/*        Form P**T, determined by a call to DGEBRD to reduce a k-by-n */
/*        matrix */

#line 300 "dorgbr.f"
	if (*k < *n) {

/*           If k < n, assume k <= m <= n */

#line 304 "dorgbr.f"
	    dorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 306 "dorgbr.f"
	} else {

/*           If k >= n, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           row downward, and set the first row and column of P**T to */
/*           those of the unit matrix */

#line 314 "dorgbr.f"
	    a[a_dim1 + 1] = 1.;
#line 315 "dorgbr.f"
	    i__1 = *n;
#line 315 "dorgbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 316 "dorgbr.f"
		a[i__ + a_dim1] = 0.;
#line 317 "dorgbr.f"
/* L40: */
#line 317 "dorgbr.f"
	    }
#line 318 "dorgbr.f"
	    i__1 = *n;
#line 318 "dorgbr.f"
	    for (j = 2; j <= i__1; ++j) {
#line 319 "dorgbr.f"
		for (i__ = j - 1; i__ >= 2; --i__) {
#line 320 "dorgbr.f"
		    a[i__ + j * a_dim1] = a[i__ - 1 + j * a_dim1];
#line 321 "dorgbr.f"
/* L50: */
#line 321 "dorgbr.f"
		}
#line 322 "dorgbr.f"
		a[j * a_dim1 + 1] = 0.;
#line 323 "dorgbr.f"
/* L60: */
#line 323 "dorgbr.f"
	    }
#line 324 "dorgbr.f"
	    if (*n > 1) {

/*              Form P**T(2:n,2:n) */

#line 328 "dorgbr.f"
		i__1 = *n - 1;
#line 328 "dorgbr.f"
		i__2 = *n - 1;
#line 328 "dorgbr.f"
		i__3 = *n - 1;
#line 328 "dorgbr.f"
		dorglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 330 "dorgbr.f"
	    }
#line 331 "dorgbr.f"
	}
#line 332 "dorgbr.f"
    }
#line 333 "dorgbr.f"
    work[1] = (doublereal) lwkopt;
#line 334 "dorgbr.f"
    return 0;

/*     End of DORGBR */

} /* dorgbr_ */

