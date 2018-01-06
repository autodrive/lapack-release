#line 1 "zgebd2.f"
/* zgebd2.f -- translated by f2c (version 20100827).
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

#line 1 "zgebd2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZGEBD2 reduces a general matrix to bidiagonal form using an unblocked algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEBD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgebd2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgebd2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgebd2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       COMPLEX*16         A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEBD2 reduces a complex general m by n matrix A to upper or lower */
/* > real bidiagonal form B by a unitary transformation: Q**H * A * P = B. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the m by n general matrix to be reduced. */
/* >          On exit, */
/* >          if m >= n, the diagonal and the first superdiagonal are */
/* >            overwritten with the upper bidiagonal matrix B; the */
/* >            elements below the diagonal, with the array TAUQ, represent */
/* >            the unitary matrix Q as a product of elementary */
/* >            reflectors, and the elements above the first superdiagonal, */
/* >            with the array TAUP, represent the unitary matrix P as */
/* >            a product of elementary reflectors; */
/* >          if m < n, the diagonal and the first subdiagonal are */
/* >            overwritten with the lower bidiagonal matrix B; the */
/* >            elements below the first subdiagonal, with the array TAUQ, */
/* >            represent the unitary matrix Q as a product of */
/* >            elementary reflectors, and the elements above the diagonal, */
/* >            with the array TAUP, represent the unitary matrix P as */
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
/* >          TAUQ is COMPLEX*16 array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the unitary matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* >          TAUP is COMPLEX*16 array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the unitary matrix P. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (max(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup complex16GEcomputational */

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
/* >     H(i) = I - tauq * v * v**H  and G(i) = I - taup * u * u**H */
/* > */
/* >  where tauq and taup are complex scalars, and v and u are complex */
/* >  vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in */
/* >  A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in */
/* >  A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* >  If m < n, */
/* > */
/* >     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m) */
/* > */
/* >  Each H(i) and G(i) has the form: */
/* > */
/* >     H(i) = I - tauq * v * v**H  and G(i) = I - taup * u * u**H */
/* > */
/* >  where tauq and taup are complex scalars, v and u are complex vectors; */
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
/* Subroutine */ int zgebd2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq, 
	doublecomplex *taup, doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublecomplex alpha;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), zlarfg_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), zlacgv_(integer *, doublecomplex *, 
	    integer *);


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

#line 226 "zgebd2.f"
    /* Parameter adjustments */
#line 226 "zgebd2.f"
    a_dim1 = *lda;
#line 226 "zgebd2.f"
    a_offset = 1 + a_dim1;
#line 226 "zgebd2.f"
    a -= a_offset;
#line 226 "zgebd2.f"
    --d__;
#line 226 "zgebd2.f"
    --e;
#line 226 "zgebd2.f"
    --tauq;
#line 226 "zgebd2.f"
    --taup;
#line 226 "zgebd2.f"
    --work;
#line 226 "zgebd2.f"

#line 226 "zgebd2.f"
    /* Function Body */
#line 226 "zgebd2.f"
    *info = 0;
#line 227 "zgebd2.f"
    if (*m < 0) {
#line 228 "zgebd2.f"
	*info = -1;
#line 229 "zgebd2.f"
    } else if (*n < 0) {
#line 230 "zgebd2.f"
	*info = -2;
#line 231 "zgebd2.f"
    } else if (*lda < max(1,*m)) {
#line 232 "zgebd2.f"
	*info = -4;
#line 233 "zgebd2.f"
    }
#line 234 "zgebd2.f"
    if (*info < 0) {
#line 235 "zgebd2.f"
	i__1 = -(*info);
#line 235 "zgebd2.f"
	xerbla_("ZGEBD2", &i__1, (ftnlen)6);
#line 236 "zgebd2.f"
	return 0;
#line 237 "zgebd2.f"
    }

#line 239 "zgebd2.f"
    if (*m >= *n) {

/*        Reduce to upper bidiagonal form */

#line 243 "zgebd2.f"
	i__1 = *n;
#line 243 "zgebd2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

#line 247 "zgebd2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 247 "zgebd2.f"
	    alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 248 "zgebd2.f"
	    i__2 = *m - i__ + 1;
/* Computing MIN */
#line 248 "zgebd2.f"
	    i__3 = i__ + 1;
#line 248 "zgebd2.f"
	    zlarfg_(&i__2, &alpha, &a[min(i__3,*m) + i__ * a_dim1], &c__1, &
		    tauq[i__]);
#line 250 "zgebd2.f"
	    i__2 = i__;
#line 250 "zgebd2.f"
	    d__[i__2] = alpha.r;
#line 251 "zgebd2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 251 "zgebd2.f"
	    a[i__2].r = 1., a[i__2].i = 0.;

/*           Apply H(i)**H to A(i:m,i+1:n) from the left */

#line 255 "zgebd2.f"
	    if (i__ < *n) {
#line 255 "zgebd2.f"
		i__2 = *m - i__ + 1;
#line 255 "zgebd2.f"
		i__3 = *n - i__;
#line 255 "zgebd2.f"
		d_cnjg(&z__1, &tauq[i__]);
#line 255 "zgebd2.f"
		zlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			z__1, &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (
			ftnlen)4);
#line 255 "zgebd2.f"
	    }
#line 258 "zgebd2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 258 "zgebd2.f"
	    i__3 = i__;
#line 258 "zgebd2.f"
	    a[i__2].r = d__[i__3], a[i__2].i = 0.;

#line 260 "zgebd2.f"
	    if (i__ < *n) {

/*              Generate elementary reflector G(i) to annihilate */
/*              A(i,i+2:n) */

#line 265 "zgebd2.f"
		i__2 = *n - i__;
#line 265 "zgebd2.f"
		zlacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 266 "zgebd2.f"
		i__2 = i__ + (i__ + 1) * a_dim1;
#line 266 "zgebd2.f"
		alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 267 "zgebd2.f"
		i__2 = *n - i__;
/* Computing MIN */
#line 267 "zgebd2.f"
		i__3 = i__ + 2;
#line 267 "zgebd2.f"
		zlarfg_(&i__2, &alpha, &a[i__ + min(i__3,*n) * a_dim1], lda, &
			taup[i__]);
#line 269 "zgebd2.f"
		i__2 = i__;
#line 269 "zgebd2.f"
		e[i__2] = alpha.r;
#line 270 "zgebd2.f"
		i__2 = i__ + (i__ + 1) * a_dim1;
#line 270 "zgebd2.f"
		a[i__2].r = 1., a[i__2].i = 0.;

/*              Apply G(i) to A(i+1:m,i+1:n) from the right */

#line 274 "zgebd2.f"
		i__2 = *m - i__;
#line 274 "zgebd2.f"
		i__3 = *n - i__;
#line 274 "zgebd2.f"
		zlarf_("Right", &i__2, &i__3, &a[i__ + (i__ + 1) * a_dim1], 
			lda, &taup[i__], &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda, &work[1], (ftnlen)5);
#line 276 "zgebd2.f"
		i__2 = *n - i__;
#line 276 "zgebd2.f"
		zlacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 277 "zgebd2.f"
		i__2 = i__ + (i__ + 1) * a_dim1;
#line 277 "zgebd2.f"
		i__3 = i__;
#line 277 "zgebd2.f"
		a[i__2].r = e[i__3], a[i__2].i = 0.;
#line 278 "zgebd2.f"
	    } else {
#line 279 "zgebd2.f"
		i__2 = i__;
#line 279 "zgebd2.f"
		taup[i__2].r = 0., taup[i__2].i = 0.;
#line 280 "zgebd2.f"
	    }
#line 281 "zgebd2.f"
/* L10: */
#line 281 "zgebd2.f"
	}
#line 282 "zgebd2.f"
    } else {

/*        Reduce to lower bidiagonal form */

#line 286 "zgebd2.f"
	i__1 = *m;
#line 286 "zgebd2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector G(i) to annihilate A(i,i+1:n) */

#line 290 "zgebd2.f"
	    i__2 = *n - i__ + 1;
#line 290 "zgebd2.f"
	    zlacgv_(&i__2, &a[i__ + i__ * a_dim1], lda);
#line 291 "zgebd2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 291 "zgebd2.f"
	    alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 292 "zgebd2.f"
	    i__2 = *n - i__ + 1;
/* Computing MIN */
#line 292 "zgebd2.f"
	    i__3 = i__ + 1;
#line 292 "zgebd2.f"
	    zlarfg_(&i__2, &alpha, &a[i__ + min(i__3,*n) * a_dim1], lda, &
		    taup[i__]);
#line 294 "zgebd2.f"
	    i__2 = i__;
#line 294 "zgebd2.f"
	    d__[i__2] = alpha.r;
#line 295 "zgebd2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 295 "zgebd2.f"
	    a[i__2].r = 1., a[i__2].i = 0.;

/*           Apply G(i) to A(i+1:m,i:n) from the right */

#line 299 "zgebd2.f"
	    if (i__ < *m) {
#line 299 "zgebd2.f"
		i__2 = *m - i__;
#line 299 "zgebd2.f"
		i__3 = *n - i__ + 1;
#line 299 "zgebd2.f"
		zlarf_("Right", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &
			taup[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1], 
			(ftnlen)5);
#line 299 "zgebd2.f"
	    }
#line 302 "zgebd2.f"
	    i__2 = *n - i__ + 1;
#line 302 "zgebd2.f"
	    zlacgv_(&i__2, &a[i__ + i__ * a_dim1], lda);
#line 303 "zgebd2.f"
	    i__2 = i__ + i__ * a_dim1;
#line 303 "zgebd2.f"
	    i__3 = i__;
#line 303 "zgebd2.f"
	    a[i__2].r = d__[i__3], a[i__2].i = 0.;

#line 305 "zgebd2.f"
	    if (i__ < *m) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(i+2:m,i) */

#line 310 "zgebd2.f"
		i__2 = i__ + 1 + i__ * a_dim1;
#line 310 "zgebd2.f"
		alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 311 "zgebd2.f"
		i__2 = *m - i__;
/* Computing MIN */
#line 311 "zgebd2.f"
		i__3 = i__ + 2;
#line 311 "zgebd2.f"
		zlarfg_(&i__2, &alpha, &a[min(i__3,*m) + i__ * a_dim1], &c__1,
			 &tauq[i__]);
#line 313 "zgebd2.f"
		i__2 = i__;
#line 313 "zgebd2.f"
		e[i__2] = alpha.r;
#line 314 "zgebd2.f"
		i__2 = i__ + 1 + i__ * a_dim1;
#line 314 "zgebd2.f"
		a[i__2].r = 1., a[i__2].i = 0.;

/*              Apply H(i)**H to A(i+1:m,i+1:n) from the left */

#line 318 "zgebd2.f"
		i__2 = *m - i__;
#line 318 "zgebd2.f"
		i__3 = *n - i__;
#line 318 "zgebd2.f"
		d_cnjg(&z__1, &tauq[i__]);
#line 318 "zgebd2.f"
		zlarf_("Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &
			c__1, &z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &
			work[1], (ftnlen)4);
#line 321 "zgebd2.f"
		i__2 = i__ + 1 + i__ * a_dim1;
#line 321 "zgebd2.f"
		i__3 = i__;
#line 321 "zgebd2.f"
		a[i__2].r = e[i__3], a[i__2].i = 0.;
#line 322 "zgebd2.f"
	    } else {
#line 323 "zgebd2.f"
		i__2 = i__;
#line 323 "zgebd2.f"
		tauq[i__2].r = 0., tauq[i__2].i = 0.;
#line 324 "zgebd2.f"
	    }
#line 325 "zgebd2.f"
/* L20: */
#line 325 "zgebd2.f"
	}
#line 326 "zgebd2.f"
    }
#line 327 "zgebd2.f"
    return 0;

/*     End of ZGEBD2 */

} /* zgebd2_ */

