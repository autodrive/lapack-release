#line 1 "dgebd2.f"
/* dgebd2.f -- translated by f2c (version 20100827).
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

#line 1 "dgebd2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DGEBD2 reduces a general matrix to bidiagonal form using an unblocked algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEBD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebd2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebd2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebd2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ), */
/*      $                   TAUQ( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEBD2 reduces a real general m by n matrix A to upper or lower */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the m by n general matrix to be reduced. */
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
/* >          D is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The diagonal elements of the bidiagonal matrix B: */
/* >          D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (min(M,N)-1) */
/* >          The off-diagonal elements of the bidiagonal matrix B: */
/* >          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; */
/* >          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ */
/* > \verbatim */
/* >          TAUQ is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the orthogonal matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* >          TAUP is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the orthogonal matrix P. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (max(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit. */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup doubleGEcomputational */

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
/* Subroutine */ int dgebd2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	taup, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), dlarfg_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *), xerbla_(char *, integer *,
	     ftnlen);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 224 "dgebd2.f"
    /* Parameter adjustments */
#line 224 "dgebd2.f"
    a_dim1 = *lda;
#line 224 "dgebd2.f"
    a_offset = 1 + a_dim1;
#line 224 "dgebd2.f"
    a -= a_offset;
#line 224 "dgebd2.f"
    --d__;
#line 224 "dgebd2.f"
    --e;
#line 224 "dgebd2.f"
    --tauq;
#line 224 "dgebd2.f"
    --taup;
#line 224 "dgebd2.f"
    --work;
#line 224 "dgebd2.f"

#line 224 "dgebd2.f"
    /* Function Body */
#line 224 "dgebd2.f"
    *info = 0;
#line 225 "dgebd2.f"
    if (*m < 0) {
#line 226 "dgebd2.f"
	*info = -1;
#line 227 "dgebd2.f"
    } else if (*n < 0) {
#line 228 "dgebd2.f"
	*info = -2;
#line 229 "dgebd2.f"
    } else if (*lda < max(1,*m)) {
#line 230 "dgebd2.f"
	*info = -4;
#line 231 "dgebd2.f"
    }
#line 232 "dgebd2.f"
    if (*info < 0) {
#line 233 "dgebd2.f"
	i__1 = -(*info);
#line 233 "dgebd2.f"
	xerbla_("DGEBD2", &i__1, (ftnlen)6);
#line 234 "dgebd2.f"
	return 0;
#line 235 "dgebd2.f"
    }

#line 237 "dgebd2.f"
    if (*m >= *n) {

/*        Reduce to upper bidiagonal form */

#line 241 "dgebd2.f"
	i__1 = *n;
#line 241 "dgebd2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

#line 245 "dgebd2.f"
	    i__2 = *m - i__ + 1;
/* Computing MIN */
#line 245 "dgebd2.f"
	    i__3 = i__ + 1;
#line 245 "dgebd2.f"
	    dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3,*m) + i__ * 
		    a_dim1], &c__1, &tauq[i__]);
#line 247 "dgebd2.f"
	    d__[i__] = a[i__ + i__ * a_dim1];
#line 248 "dgebd2.f"
	    a[i__ + i__ * a_dim1] = 1.;

/*           Apply H(i) to A(i:m,i+1:n) from the left */

#line 252 "dgebd2.f"
	    if (i__ < *n) {
#line 252 "dgebd2.f"
		i__2 = *m - i__ + 1;
#line 252 "dgebd2.f"
		i__3 = *n - i__;
#line 252 "dgebd2.f"
		dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			tauq[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]
			, (ftnlen)4);
#line 252 "dgebd2.f"
	    }
#line 255 "dgebd2.f"
	    a[i__ + i__ * a_dim1] = d__[i__];

#line 257 "dgebd2.f"
	    if (i__ < *n) {

/*              Generate elementary reflector G(i) to annihilate */
/*              A(i,i+2:n) */

#line 262 "dgebd2.f"
		i__2 = *n - i__;
/* Computing MIN */
#line 262 "dgebd2.f"
		i__3 = i__ + 2;
#line 262 "dgebd2.f"
		dlarfg_(&i__2, &a[i__ + (i__ + 1) * a_dim1], &a[i__ + min(
			i__3,*n) * a_dim1], lda, &taup[i__]);
#line 264 "dgebd2.f"
		e[i__] = a[i__ + (i__ + 1) * a_dim1];
#line 265 "dgebd2.f"
		a[i__ + (i__ + 1) * a_dim1] = 1.;

/*              Apply G(i) to A(i+1:m,i+1:n) from the right */

#line 269 "dgebd2.f"
		i__2 = *m - i__;
#line 269 "dgebd2.f"
		i__3 = *n - i__;
#line 269 "dgebd2.f"
		dlarf_("Right", &i__2, &i__3, &a[i__ + (i__ + 1) * a_dim1], 
			lda, &taup[i__], &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda, &work[1], (ftnlen)5);
#line 271 "dgebd2.f"
		a[i__ + (i__ + 1) * a_dim1] = e[i__];
#line 272 "dgebd2.f"
	    } else {
#line 273 "dgebd2.f"
		taup[i__] = 0.;
#line 274 "dgebd2.f"
	    }
#line 275 "dgebd2.f"
/* L10: */
#line 275 "dgebd2.f"
	}
#line 276 "dgebd2.f"
    } else {

/*        Reduce to lower bidiagonal form */

#line 280 "dgebd2.f"
	i__1 = *m;
#line 280 "dgebd2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector G(i) to annihilate A(i,i+1:n) */

#line 284 "dgebd2.f"
	    i__2 = *n - i__ + 1;
/* Computing MIN */
#line 284 "dgebd2.f"
	    i__3 = i__ + 1;
#line 284 "dgebd2.f"
	    dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + min(i__3,*n) * 
		    a_dim1], lda, &taup[i__]);
#line 286 "dgebd2.f"
	    d__[i__] = a[i__ + i__ * a_dim1];
#line 287 "dgebd2.f"
	    a[i__ + i__ * a_dim1] = 1.;

/*           Apply G(i) to A(i+1:m,i:n) from the right */

#line 291 "dgebd2.f"
	    if (i__ < *m) {
#line 291 "dgebd2.f"
		i__2 = *m - i__;
#line 291 "dgebd2.f"
		i__3 = *n - i__ + 1;
#line 291 "dgebd2.f"
		dlarf_("Right", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &
			taup[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1], 
			(ftnlen)5);
#line 291 "dgebd2.f"
	    }
#line 294 "dgebd2.f"
	    a[i__ + i__ * a_dim1] = d__[i__];

#line 296 "dgebd2.f"
	    if (i__ < *m) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(i+2:m,i) */

#line 301 "dgebd2.f"
		i__2 = *m - i__;
/* Computing MIN */
#line 301 "dgebd2.f"
		i__3 = i__ + 2;
#line 301 "dgebd2.f"
		dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*m) + 
			i__ * a_dim1], &c__1, &tauq[i__]);
#line 303 "dgebd2.f"
		e[i__] = a[i__ + 1 + i__ * a_dim1];
#line 304 "dgebd2.f"
		a[i__ + 1 + i__ * a_dim1] = 1.;

/*              Apply H(i) to A(i+1:m,i+1:n) from the left */

#line 308 "dgebd2.f"
		i__2 = *m - i__;
#line 308 "dgebd2.f"
		i__3 = *n - i__;
#line 308 "dgebd2.f"
		dlarf_("Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &
			c__1, &tauq[i__], &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda, &work[1], (ftnlen)4);
#line 310 "dgebd2.f"
		a[i__ + 1 + i__ * a_dim1] = e[i__];
#line 311 "dgebd2.f"
	    } else {
#line 312 "dgebd2.f"
		tauq[i__] = 0.;
#line 313 "dgebd2.f"
	    }
#line 314 "dgebd2.f"
/* L20: */
#line 314 "dgebd2.f"
	}
#line 315 "dgebd2.f"
    }
#line 316 "dgebd2.f"
    return 0;

/*     End of DGEBD2 */

} /* dgebd2_ */

