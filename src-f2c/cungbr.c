#line 1 "cungbr.f"
/* cungbr.f -- translated by f2c (version 20100827).
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

#line 1 "cungbr.f"
/* Table of constant values */

static integer c_n1 = -1;

/* > \brief \b CUNGBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNGBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cungbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cungbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cungbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          VECT */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNGBR generates one of the complex unitary matrices Q or P**H */
/* > determined by CGEBRD when reducing a complex matrix A to bidiagonal */
/* > form: A = Q * B * P**H.  Q and P**H are defined as products of */
/* > elementary reflectors H(i) or G(i) respectively. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q */
/* > is of order M: */
/* > if m >= k, Q = H(1) H(2) . . . H(k) and CUNGBR returns the first n */
/* > columns of Q, where m >= n >= k; */
/* > if m < k, Q = H(1) H(2) . . . H(m-1) and CUNGBR returns Q as an */
/* > M-by-M matrix. */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**H */
/* > is of order N: */
/* > if k < n, P**H = G(k) . . . G(2) G(1) and CUNGBR returns the first m */
/* > rows of P**H, where n >= m >= k; */
/* > if k >= n, P**H = G(n-1) . . . G(2) G(1) and CUNGBR returns P**H as */
/* > an N-by-N matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          Specifies whether the matrix Q or the matrix P**H is */
/* >          required, as defined in the transformation applied by CGEBRD: */
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
/* >          matrix reduced by CGEBRD. */
/* >          If VECT = 'P', the number of rows in the original K-by-N */
/* >          matrix reduced by CGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by CGEBRD. */
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
/* >          TAU is COMPLEX array, dimension */
/* >                                (min(M,K)) if VECT = 'Q' */
/* >                                (min(N,K)) if VECT = 'P' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i), which determines Q or P**H, as */
/* >          returned by CGEBRD in its array argument TAUQ or TAUP. */
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

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int cungbr_(char *vect, integer *m, integer *n, integer *k, 
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cunglq_(
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), cungqr_(
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);
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

#line 199 "cungbr.f"
    /* Parameter adjustments */
#line 199 "cungbr.f"
    a_dim1 = *lda;
#line 199 "cungbr.f"
    a_offset = 1 + a_dim1;
#line 199 "cungbr.f"
    a -= a_offset;
#line 199 "cungbr.f"
    --tau;
#line 199 "cungbr.f"
    --work;
#line 199 "cungbr.f"

#line 199 "cungbr.f"
    /* Function Body */
#line 199 "cungbr.f"
    *info = 0;
#line 200 "cungbr.f"
    wantq = lsame_(vect, "Q", (ftnlen)1, (ftnlen)1);
#line 201 "cungbr.f"
    mn = min(*m,*n);
#line 202 "cungbr.f"
    lquery = *lwork == -1;
#line 203 "cungbr.f"
    if (! wantq && ! lsame_(vect, "P", (ftnlen)1, (ftnlen)1)) {
#line 204 "cungbr.f"
	*info = -1;
#line 205 "cungbr.f"
    } else if (*m < 0) {
#line 206 "cungbr.f"
	*info = -2;
#line 207 "cungbr.f"
    } else if (*n < 0 || wantq && (*n > *m || *n < min(*m,*k)) || ! wantq && (
	    *m > *n || *m < min(*n,*k))) {
#line 210 "cungbr.f"
	*info = -3;
#line 211 "cungbr.f"
    } else if (*k < 0) {
#line 212 "cungbr.f"
	*info = -4;
#line 213 "cungbr.f"
    } else if (*lda < max(1,*m)) {
#line 214 "cungbr.f"
	*info = -6;
#line 215 "cungbr.f"
    } else if (*lwork < max(1,mn) && ! lquery) {
#line 216 "cungbr.f"
	*info = -9;
#line 217 "cungbr.f"
    }

#line 219 "cungbr.f"
    if (*info == 0) {
#line 220 "cungbr.f"
	work[1].r = 1., work[1].i = 0.;
#line 221 "cungbr.f"
	if (wantq) {
#line 222 "cungbr.f"
	    if (*m >= *k) {
#line 223 "cungbr.f"
		cungqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 224 "cungbr.f"
	    } else {
#line 225 "cungbr.f"
		if (*m > 1) {
#line 226 "cungbr.f"
		    i__1 = *m - 1;
#line 226 "cungbr.f"
		    i__2 = *m - 1;
#line 226 "cungbr.f"
		    i__3 = *m - 1;
#line 226 "cungbr.f"
		    cungqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 228 "cungbr.f"
		}
#line 229 "cungbr.f"
	    }
#line 230 "cungbr.f"
	} else {
#line 231 "cungbr.f"
	    if (*k < *n) {
#line 232 "cungbr.f"
		cunglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
#line 233 "cungbr.f"
	    } else {
#line 234 "cungbr.f"
		if (*n > 1) {
#line 235 "cungbr.f"
		    i__1 = *n - 1;
#line 235 "cungbr.f"
		    i__2 = *n - 1;
#line 235 "cungbr.f"
		    i__3 = *n - 1;
#line 235 "cungbr.f"
		    cunglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &
			    tau[1], &work[1], &c_n1, &iinfo);
#line 237 "cungbr.f"
		}
#line 238 "cungbr.f"
	    }
#line 239 "cungbr.f"
	}
#line 240 "cungbr.f"
	lwkopt = (integer) work[1].r;
#line 241 "cungbr.f"
	lwkopt = max(lwkopt,mn);
#line 242 "cungbr.f"
    }

#line 244 "cungbr.f"
    if (*info != 0) {
#line 245 "cungbr.f"
	i__1 = -(*info);
#line 245 "cungbr.f"
	xerbla_("CUNGBR", &i__1, (ftnlen)6);
#line 246 "cungbr.f"
	return 0;
#line 247 "cungbr.f"
    } else if (lquery) {
#line 248 "cungbr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 249 "cungbr.f"
	return 0;
#line 250 "cungbr.f"
    }

/*     Quick return if possible */

#line 254 "cungbr.f"
    if (*m == 0 || *n == 0) {
#line 255 "cungbr.f"
	work[1].r = 1., work[1].i = 0.;
#line 256 "cungbr.f"
	return 0;
#line 257 "cungbr.f"
    }

#line 259 "cungbr.f"
    if (wantq) {

/*        Form Q, determined by a call to CGEBRD to reduce an m-by-k */
/*        matrix */

#line 264 "cungbr.f"
	if (*m >= *k) {

/*           If m >= k, assume m >= n >= k */

#line 268 "cungbr.f"
	    cungqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 270 "cungbr.f"
	} else {

/*           If m < k, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           column to the right, and set the first row and column of Q */
/*           to those of the unit matrix */

#line 278 "cungbr.f"
	    for (j = *m; j >= 2; --j) {
#line 279 "cungbr.f"
		i__1 = j * a_dim1 + 1;
#line 279 "cungbr.f"
		a[i__1].r = 0., a[i__1].i = 0.;
#line 280 "cungbr.f"
		i__1 = *m;
#line 280 "cungbr.f"
		for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 281 "cungbr.f"
		    i__2 = i__ + j * a_dim1;
#line 281 "cungbr.f"
		    i__3 = i__ + (j - 1) * a_dim1;
#line 281 "cungbr.f"
		    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 282 "cungbr.f"
/* L10: */
#line 282 "cungbr.f"
		}
#line 283 "cungbr.f"
/* L20: */
#line 283 "cungbr.f"
	    }
#line 284 "cungbr.f"
	    i__1 = a_dim1 + 1;
#line 284 "cungbr.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 285 "cungbr.f"
	    i__1 = *m;
#line 285 "cungbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 286 "cungbr.f"
		i__2 = i__ + a_dim1;
#line 286 "cungbr.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 287 "cungbr.f"
/* L30: */
#line 287 "cungbr.f"
	    }
#line 288 "cungbr.f"
	    if (*m > 1) {

/*              Form Q(2:m,2:m) */

#line 292 "cungbr.f"
		i__1 = *m - 1;
#line 292 "cungbr.f"
		i__2 = *m - 1;
#line 292 "cungbr.f"
		i__3 = *m - 1;
#line 292 "cungbr.f"
		cungqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 294 "cungbr.f"
	    }
#line 295 "cungbr.f"
	}
#line 296 "cungbr.f"
    } else {

/*        Form P**H, determined by a call to CGEBRD to reduce a k-by-n */
/*        matrix */

#line 301 "cungbr.f"
	if (*k < *n) {

/*           If k < n, assume k <= m <= n */

#line 305 "cungbr.f"
	    cunglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

#line 307 "cungbr.f"
	} else {

/*           If k >= n, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           row downward, and set the first row and column of P**H to */
/*           those of the unit matrix */

#line 315 "cungbr.f"
	    i__1 = a_dim1 + 1;
#line 315 "cungbr.f"
	    a[i__1].r = 1., a[i__1].i = 0.;
#line 316 "cungbr.f"
	    i__1 = *n;
#line 316 "cungbr.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 317 "cungbr.f"
		i__2 = i__ + a_dim1;
#line 317 "cungbr.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 318 "cungbr.f"
/* L40: */
#line 318 "cungbr.f"
	    }
#line 319 "cungbr.f"
	    i__1 = *n;
#line 319 "cungbr.f"
	    for (j = 2; j <= i__1; ++j) {
#line 320 "cungbr.f"
		for (i__ = j - 1; i__ >= 2; --i__) {
#line 321 "cungbr.f"
		    i__2 = i__ + j * a_dim1;
#line 321 "cungbr.f"
		    i__3 = i__ - 1 + j * a_dim1;
#line 321 "cungbr.f"
		    a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
#line 322 "cungbr.f"
/* L50: */
#line 322 "cungbr.f"
		}
#line 323 "cungbr.f"
		i__2 = j * a_dim1 + 1;
#line 323 "cungbr.f"
		a[i__2].r = 0., a[i__2].i = 0.;
#line 324 "cungbr.f"
/* L60: */
#line 324 "cungbr.f"
	    }
#line 325 "cungbr.f"
	    if (*n > 1) {

/*              Form P**H(2:n,2:n) */

#line 329 "cungbr.f"
		i__1 = *n - 1;
#line 329 "cungbr.f"
		i__2 = *n - 1;
#line 329 "cungbr.f"
		i__3 = *n - 1;
#line 329 "cungbr.f"
		cunglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
#line 331 "cungbr.f"
	    }
#line 332 "cungbr.f"
	}
#line 333 "cungbr.f"
    }
#line 334 "cungbr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 335 "cungbr.f"
    return 0;

/*     End of CUNGBR */

} /* cungbr_ */

