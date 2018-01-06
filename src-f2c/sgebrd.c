#line 1 "sgebrd.f"
/* sgebrd.f -- translated by f2c (version 20100827).
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

#line 1 "sgebrd.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static doublereal c_b21 = -1.;
static doublereal c_b22 = 1.;

/* > \brief \b SGEBRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEBRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgebrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgebrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgebrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), D( * ), E( * ), TAUP( * ), */
/*      $                   TAUQ( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEBRD reduces a general real M-by-N matrix A to upper or lower */
/* > bidiagonal form B by an orthogonal transformation: Q**T * A * P = B. */
/* > */
/* > If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows in the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns in the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M-by-N general matrix to be reduced. */
/* >          On exit, */
/* >          if m >= n, the diagonal and the first superdiagonal are */
/* >            overwritten with the upper bidiagonal matrix B; the */
/* >            elements below the diagonal, with the array TAUQ, represent */
/* >            the orthogonal matrix Q as a product of elementary */
/* >            reflectors, and the elements above the first superdiagonal, */
/* >            with the array TAUP, represent the orthogonal matrix P as */
/* >            a product of elementary reflectors; */
/* >          if m < n, the diagonal and the first subdiagonal are */
/* >            overwritten with the lower bidiagonal matrix B; the */
/* >            elements below the first subdiagonal, with the array TAUQ, */
/* >            represent the orthogonal matrix Q as a product of */
/* >            elementary reflectors, and the elements above the diagonal, */
/* >            with the array TAUP, represent the orthogonal matrix P as */
/* >            a product of elementary reflectors. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (min(M,N)) */
/* >          The diagonal elements of the bidiagonal matrix B: */
/* >          D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (min(M,N)-1) */
/* >          The off-diagonal elements of the bidiagonal matrix B: */
/* >          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; */
/* >          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ */
/* > \verbatim */
/* >          TAUQ is REAL array dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the orthogonal matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* >          TAUP is REAL array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the orthogonal matrix P. See Further Details. */
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
/* >          The length of the array WORK.  LWORK >= max(1,M,N). */
/* >          For optimum performance LWORK >= (M+N)*NB, where NB */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrices Q and P are represented as products of elementary */
/* >  reflectors: */
/* > */
/* >  If m >= n, */
/* > */
/* >     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1) */
/* > */
/* >  Each H(i) and G(i) has the form: */
/* > */
/* >     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T */
/* > */
/* >  where tauq and taup are real scalars, and v and u are real vectors; */
/* >  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i); */
/* >  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n); */
/* >  tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* >  If m < n, */
/* > */
/* >     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m) */
/* > */
/* >  Each H(i) and G(i) has the form: */
/* > */
/* >     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T */
/* > */
/* >  where tauq and taup are real scalars, and v and u are real vectors; */
/* >  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i); */
/* >  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n); */
/* >  tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* >  The contents of A on exit are illustrated by the following examples: */
/* > */
/* >  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n): */
/* > */
/* >    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 ) */
/* >    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 ) */
/* >    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 ) */
/* >    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 ) */
/* >    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 ) */
/* >    (  v1  v2  v3  v4  v5 ) */
/* > */
/* >  where d and e denote diagonal and off-diagonal elements of B, vi */
/* >  denotes an element of the vector defining H(i), and ui an element of */
/* >  the vector defining G(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgebrd_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	taup, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, nb, nx;
    static doublereal ws;
    static integer nbmin, iinfo;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer minmn;
    extern /* Subroutine */ int sgebd2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *), slabrd_(integer *, integer *, integer *
	    , doublereal *, integer *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwrkx, ldwrky, lwkopt;
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

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 247 "sgebrd.f"
    /* Parameter adjustments */
#line 247 "sgebrd.f"
    a_dim1 = *lda;
#line 247 "sgebrd.f"
    a_offset = 1 + a_dim1;
#line 247 "sgebrd.f"
    a -= a_offset;
#line 247 "sgebrd.f"
    --d__;
#line 247 "sgebrd.f"
    --e;
#line 247 "sgebrd.f"
    --tauq;
#line 247 "sgebrd.f"
    --taup;
#line 247 "sgebrd.f"
    --work;
#line 247 "sgebrd.f"

#line 247 "sgebrd.f"
    /* Function Body */
#line 247 "sgebrd.f"
    *info = 0;
/* Computing MAX */
#line 248 "sgebrd.f"
    i__1 = 1, i__2 = ilaenv_(&c__1, "SGEBRD", " ", m, n, &c_n1, &c_n1, (
	    ftnlen)6, (ftnlen)1);
#line 248 "sgebrd.f"
    nb = max(i__1,i__2);
#line 249 "sgebrd.f"
    lwkopt = (*m + *n) * nb;
#line 250 "sgebrd.f"
    work[1] = (doublereal) lwkopt;
#line 251 "sgebrd.f"
    lquery = *lwork == -1;
#line 252 "sgebrd.f"
    if (*m < 0) {
#line 253 "sgebrd.f"
	*info = -1;
#line 254 "sgebrd.f"
    } else if (*n < 0) {
#line 255 "sgebrd.f"
	*info = -2;
#line 256 "sgebrd.f"
    } else if (*lda < max(1,*m)) {
#line 257 "sgebrd.f"
	*info = -4;
#line 258 "sgebrd.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 258 "sgebrd.f"
	i__1 = max(1,*m);
#line 258 "sgebrd.f"
	if (*lwork < max(i__1,*n) && ! lquery) {
#line 259 "sgebrd.f"
	    *info = -10;
#line 260 "sgebrd.f"
	}
#line 260 "sgebrd.f"
    }
#line 261 "sgebrd.f"
    if (*info < 0) {
#line 262 "sgebrd.f"
	i__1 = -(*info);
#line 262 "sgebrd.f"
	xerbla_("SGEBRD", &i__1, (ftnlen)6);
#line 263 "sgebrd.f"
	return 0;
#line 264 "sgebrd.f"
    } else if (lquery) {
#line 265 "sgebrd.f"
	return 0;
#line 266 "sgebrd.f"
    }

/*     Quick return if possible */

#line 270 "sgebrd.f"
    minmn = min(*m,*n);
#line 271 "sgebrd.f"
    if (minmn == 0) {
#line 272 "sgebrd.f"
	work[1] = 1.;
#line 273 "sgebrd.f"
	return 0;
#line 274 "sgebrd.f"
    }

#line 276 "sgebrd.f"
    ws = (doublereal) max(*m,*n);
#line 277 "sgebrd.f"
    ldwrkx = *m;
#line 278 "sgebrd.f"
    ldwrky = *n;

#line 280 "sgebrd.f"
    if (nb > 1 && nb < minmn) {

/*        Set the crossover point NX. */

/* Computing MAX */
#line 284 "sgebrd.f"
	i__1 = nb, i__2 = ilaenv_(&c__3, "SGEBRD", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 284 "sgebrd.f"
	nx = max(i__1,i__2);

/*        Determine when to switch from blocked to unblocked code. */

#line 288 "sgebrd.f"
	if (nx < minmn) {
#line 289 "sgebrd.f"
	    ws = (doublereal) ((*m + *n) * nb);
#line 290 "sgebrd.f"
	    if ((doublereal) (*lwork) < ws) {

/*              Not enough work space for the optimal NB, consider using */
/*              a smaller block size. */

#line 295 "sgebrd.f"
		nbmin = ilaenv_(&c__2, "SGEBRD", " ", m, n, &c_n1, &c_n1, (
			ftnlen)6, (ftnlen)1);
#line 296 "sgebrd.f"
		if (*lwork >= (*m + *n) * nbmin) {
#line 297 "sgebrd.f"
		    nb = *lwork / (*m + *n);
#line 298 "sgebrd.f"
		} else {
#line 299 "sgebrd.f"
		    nb = 1;
#line 300 "sgebrd.f"
		    nx = minmn;
#line 301 "sgebrd.f"
		}
#line 302 "sgebrd.f"
	    }
#line 303 "sgebrd.f"
	}
#line 304 "sgebrd.f"
    } else {
#line 305 "sgebrd.f"
	nx = minmn;
#line 306 "sgebrd.f"
    }

#line 308 "sgebrd.f"
    i__1 = minmn - nx;
#line 308 "sgebrd.f"
    i__2 = nb;
#line 308 "sgebrd.f"
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*        Reduce rows and columns i:i+nb-1 to bidiagonal form and return */
/*        the matrices X and Y which are needed to update the unreduced */
/*        part of the matrix */

#line 314 "sgebrd.f"
	i__3 = *m - i__ + 1;
#line 314 "sgebrd.f"
	i__4 = *n - i__ + 1;
#line 314 "sgebrd.f"
	slabrd_(&i__3, &i__4, &nb, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[
		i__], &tauq[i__], &taup[i__], &work[1], &ldwrkx, &work[ldwrkx 
		* nb + 1], &ldwrky);

/*        Update the trailing submatrix A(i+nb:m,i+nb:n), using an update */
/*        of the form  A := A - V*Y**T - X*U**T */

#line 321 "sgebrd.f"
	i__3 = *m - i__ - nb + 1;
#line 321 "sgebrd.f"
	i__4 = *n - i__ - nb + 1;
#line 321 "sgebrd.f"
	sgemm_("No transpose", "Transpose", &i__3, &i__4, &nb, &c_b21, &a[i__ 
		+ nb + i__ * a_dim1], lda, &work[ldwrkx * nb + nb + 1], &
		ldwrky, &c_b22, &a[i__ + nb + (i__ + nb) * a_dim1], lda, (
		ftnlen)12, (ftnlen)9);
#line 325 "sgebrd.f"
	i__3 = *m - i__ - nb + 1;
#line 325 "sgebrd.f"
	i__4 = *n - i__ - nb + 1;
#line 325 "sgebrd.f"
	sgemm_("No transpose", "No transpose", &i__3, &i__4, &nb, &c_b21, &
		work[nb + 1], &ldwrkx, &a[i__ + (i__ + nb) * a_dim1], lda, &
		c_b22, &a[i__ + nb + (i__ + nb) * a_dim1], lda, (ftnlen)12, (
		ftnlen)12);

/*        Copy diagonal and off-diagonal elements of B back into A */

#line 331 "sgebrd.f"
	if (*m >= *n) {
#line 332 "sgebrd.f"
	    i__3 = i__ + nb - 1;
#line 332 "sgebrd.f"
	    for (j = i__; j <= i__3; ++j) {
#line 333 "sgebrd.f"
		a[j + j * a_dim1] = d__[j];
#line 334 "sgebrd.f"
		a[j + (j + 1) * a_dim1] = e[j];
#line 335 "sgebrd.f"
/* L10: */
#line 335 "sgebrd.f"
	    }
#line 336 "sgebrd.f"
	} else {
#line 337 "sgebrd.f"
	    i__3 = i__ + nb - 1;
#line 337 "sgebrd.f"
	    for (j = i__; j <= i__3; ++j) {
#line 338 "sgebrd.f"
		a[j + j * a_dim1] = d__[j];
#line 339 "sgebrd.f"
		a[j + 1 + j * a_dim1] = e[j];
#line 340 "sgebrd.f"
/* L20: */
#line 340 "sgebrd.f"
	    }
#line 341 "sgebrd.f"
	}
#line 342 "sgebrd.f"
/* L30: */
#line 342 "sgebrd.f"
    }

/*     Use unblocked code to reduce the remainder of the matrix */

#line 346 "sgebrd.f"
    i__2 = *m - i__ + 1;
#line 346 "sgebrd.f"
    i__1 = *n - i__ + 1;
#line 346 "sgebrd.f"
    sgebd2_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], &
	    tauq[i__], &taup[i__], &work[1], &iinfo);
#line 348 "sgebrd.f"
    work[1] = ws;
#line 349 "sgebrd.f"
    return 0;

/*     End of SGEBRD */

} /* sgebrd_ */

