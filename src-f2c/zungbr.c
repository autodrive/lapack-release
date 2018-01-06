#line 1 "zungbr.f"
/* zungbr.f -- translated by f2c (version 20100827).
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

#line 1 "zungbr.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief \b ZUNGBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNGBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          VECT */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNGBR generates one of the complex unitary matrices Q or P**H */
/* > determined by ZGEBRD when reducing a complex matrix A to bidiagonal */
/* > form: A = Q * B * P**H.  Q and P**H are defined as products of */
/* > elementary reflectors H(i) or G(i) respectively. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q */
/* > is of order M: */
/* > if m >= k, Q = H(1) H(2) . . . H(k) and ZUNGBR returns the first n */
/* > columns of Q, where m >= n >= k; */
/* > if m < k, Q = H(1) H(2) . . . H(m-1) and ZUNGBR returns Q as an */
/* > M-by-M matrix. */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**H */
/* > is of order N: */
/* > if k < n, P**H = G(k) . . . G(2) G(1) and ZUNGBR returns the first m */
/* > rows of P**H, where n >= m >= k; */
/* > if k >= n, P**H = G(n-1) . . . G(2) G(1) and ZUNGBR returns P**H as */
/* > an N-by-N matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          Specifies whether the matrix Q or the matrix P**H is */
/* >          required, as defined in the transformation applied by ZGEBRD: */
/* >          = 'Q':  generate Q; */
/* >          = 'P':  generate P**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix Q or P**H to be returned. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix Q or P**H to be returned. */
/* >          N >= 0. */
/* >          If VECT = 'Q', M >= N >= min(M,K); */
/* >          if VECT = 'P', N >= M >= min(N,K). */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          If VECT = 'Q', the number of columns in the original M-by-K */
/* >          matrix reduced by ZGEBRD. */
/* >          If VECT = 'P', the number of rows in the original K-by-N */
/* >          matrix reduced by ZGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by ZGEBRD. */
/* >          On exit, the M-by-N matrix Q or P**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= M. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension */
/* >                                (min(M,K)) if VECT = 'Q' */
/* >                                (min(N,K)) if VECT = 'P' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i), which determines Q or P**H, as */
/* >          returned by ZGEBRD in its array argument TAUQ or TAUP. */
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

/* > \ingroup complex16GBcomputational */

/*  ===================================================================== */
/* Subroutine */ int zungbr_(char *vect, integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info, ftnlen vect_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, mn;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical wantq;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zunglq_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);


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

#line 199 "zungbr.f"
    /* Parameter adjustments */
#line 199 "zungbr.f"
    a_dim1 = *lda;
#line 199 "zungbr.f"
    a_offset = 1 + a_dim1;
#line 199 "zungbr.f"
    a -= a_offset;
#line 199 "zungbr.f"
    --tau;
#line 199 "zungbr.f"
    --work;
#line 199 "zungbr.f"

#line 199 "zungbr.f"
    /* Function Body */
#line 199 "zungbr.f"
    *info = 0;
#line 200 "zungbr.f"
    wantq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 201 "zungbr.f"
    mn = min(*m,*n);
#line 202 "zungbr.f"
    lquery = *lwork == -1;
#line 203 "zungbr.f"
    if (! wantq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 204 "zungbr.f"
	*info = -1;
#line 205 "zungbr.f"
    } else if (*m < 0) {
#line 206 "zungbr.f"
	*info = -2;
#line 207 "zungbr.f"
    } else if (*n < 0 || wantq && (*n > *m || *n < min(*m,*k)) || ! wantq && (
	    *m > *n || *m < min(*n,*k))) {
#line 210 "zungbr.f"
	*info = -3;
#line 211 "zungbr.f"
    } else if (*k < 0) {
#line 212 "zungbr.f"
	*info = -4;
#line 213 "zungbr.f"
    } else if (*lda < max(1,*m)) {
#line 214 "zungbr.f"
	*info = -6;
#line 215 "zungbr.f"
    } else if (*lwork < max(1,mn) && ! lquery) {
#line 216 "zungbr.f"
	*info = -9;
#line 217 "zungbr.f"
    }

#line 219 "zungbr.f"
    if (*info == 0) {
#line 220 "zungbr.f"
	work[1].r = 1., work[1].i = 0.;
#line 221 "zungbr.f"
	if (wantq) {
#line 222 "zungbr.f"
	    if (*m >= *k) {
#line 223 "zungbr.f"
		zungqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 224 "zungbr.f"
	    } else {
#line 225 "zungbr.f"
		if (*m > 1) {
#line 226 "zungbr.f"
		    i__1 = *m - 1;
#line 226 "zungbr.f"
		    i__2 = *m - 1;
#line 226 "zungbr.f"
		    i__3 = *m - 1;
#line 226 "zungbr.f"
		    zungqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 228 "zungbr.f"
		}
#line 229 "zungbr.f"
	    }
#line 230 "zungbr.f"
	} else {
#line 231 "zungbr.f"
	    if (*k < *n) {
#line 232 "zungbr.f"
		zunglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 233 "zungbr.f"
	    } else {
#line 234 "zungbr.f"
		if (*n > 1) {
#line 235 "zungbr.f"
		    i__1 = *n - 1;
#line 235 "zungbr.f"
		    i__2 = *n - 1;
#line 235 "zungbr.f"
		    i__3 = *n - 1;
#line 235 "zungbr.f"
		    zunglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 237 "zungbr.f"
		}
#line 238 "zungbr.f"
	    }
#line 239 "zungbr.f"
	}
#line 240 "zungbr.f"
	lwkopt = (integer) work[1].r;
#line 241 "zungbr.f"
	lwkopt = max(lwkopt,mn);
#line 242 "zungbr.f"
    }

#line 244 "zungbr.f"
    if (*info != 0) {
#line 245 "zungbr.f"
	i__1 = -(*info);
#line 245 "zungbr.f"
	xerbla_("ZUNGBR", &i__1, (ftnlen)6);
#line 246 "zungbr.f"
	return 0;
#line 247 "zungbr.f"
    } else if (lquery) {
#line 248 "zungbr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 249 "zungbr.f"
	return 0;
#line 250 "zungbr.f"
    }

/*     Quick return if possible */

#line 254 "zungbr.f"
    if (*m == 0 || *n == 0) {
#line 255 "zungbr.f"
	work[1].r = 1., work[1].i = 0.;
#line 256 "zungbr.f"
	return 0;
#line 257 "zungbr.f"
    }

#line 259 "zungbr.f"
    if (wantq) {

/*        Form Q, determined by a call to ZGEBRD to reduce an m-by-k */
/*        matrix */

#line 264 "zungbr.f"
	if (*m >= *k) {

/*           If m >= k, assume m >= n >= k */

#line 268 "zungbr.f"
	    zungqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 270 "zungbr.f"
	} else {

/*           If m < k, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           column to the right, and set the first row and column of Q */
/*           to those of the unit matrix */

#line 278 "zungbr.f"
	    for (j = *m; j >= 2; --j) {
#line 279 "zungbr.f"
		i__1 = j * a_dim1 + 1;
#line 279 "zungbr.f"
		a[i__1].r = 0., a[i__1].i = 0.;
#line 280 "zungbr.f"
		i__1 = *m;
#line 280 "zungbr.f"
		for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 281 "zungbr.f"
		    i__2 = i__ + j * a_dim1;
#line 281 "zungbr.f"
		    i__3 = i__ + (j - 1) * a_dim1;
#line 281 "zungbr.f"
		    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 282 "zungbr.f"
/* L10: */
#line 282 "zungbr.f"
		}
#line 283 "zungbr.f"
/* L20: */
#line 283 "zungbr.f"
	    }
#line 284 "zungbr.f"
	    i__1 = a_dim1 + 1;
#line 284 "zungbr.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 285 "zungbr.f"
	    i__1 = *m;
#line 285 "zungbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 286 "zungbr.f"
		i__2 = i__ + a_dim1;
#line 286 "zungbr.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 287 "zungbr.f"
/* L30: */
#line 287 "zungbr.f"
	    }
#line 288 "zungbr.f"
	    if (*m > 1) {

/*              Form Q(2:m,2:m) */

#line 292 "zungbr.f"
		i__1 = *m - 1;
#line 292 "zungbr.f"
		i__2 = *m - 1;
#line 292 "zungbr.f"
		i__3 = *m - 1;
#line 292 "zungbr.f"
		zungqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 294 "zungbr.f"
	    }
#line 295 "zungbr.f"
	}
#line 296 "zungbr.f"
    } else {

/*        Form P**H, determined by a call to ZGEBRD to reduce a k-by-n */
/*        matrix */

#line 301 "zungbr.f"
	if (*k < *n) {

/*           If k < n, assume k <= m <= n */

#line 305 "zungbr.f"
	    zunglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 307 "zungbr.f"
	} else {

/*           If k >= n, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           row downward, and set the first row and column of P**H to */
/*           those of the unit matrix */

#line 315 "zungbr.f"
	    i__1 = a_dim1 + 1;
#line 315 "zungbr.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 316 "zungbr.f"
	    i__1 = *n;
#line 316 "zungbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 317 "zungbr.f"
		i__2 = i__ + a_dim1;
#line 317 "zungbr.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 318 "zungbr.f"
/* L40: */
#line 318 "zungbr.f"
	    }
#line 319 "zungbr.f"
	    i__1 = *n;
#line 319 "zungbr.f"
	    for (j = 2; j <= i__1; ++j) {
#line 320 "zungbr.f"
		for (i__ = j - 1; i__ >= 2; --i__) {
#line 321 "zungbr.f"
		    i__2 = i__ + j * a_dim1;
#line 321 "zungbr.f"
		    i__3 = i__ - 1 + j * a_dim1;
#line 321 "zungbr.f"
		    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 322 "zungbr.f"
/* L50: */
#line 322 "zungbr.f"
		}
#line 323 "zungbr.f"
		i__2 = j * a_dim1 + 1;
#line 323 "zungbr.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 324 "zungbr.f"
/* L60: */
#line 324 "zungbr.f"
	    }
#line 325 "zungbr.f"
	    if (*n > 1) {

/*              Form P**H(2:n,2:n) */

#line 329 "zungbr.f"
		i__1 = *n - 1;
#line 329 "zungbr.f"
		i__2 = *n - 1;
#line 329 "zungbr.f"
		i__3 = *n - 1;
#line 329 "zungbr.f"
		zunglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 331 "zungbr.f"
	    }
#line 332 "zungbr.f"
	}
#line 333 "zungbr.f"
    }
#line 334 "zungbr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 335 "zungbr.f"
    return 0;

/*     End of ZUNGBR */

} /* zungbr_ */

