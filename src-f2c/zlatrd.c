#line 1 "zlatrd.f"
/* zlatrd.f -- translated by f2c (version 20100827).
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

#line 1 "zlatrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLATRD reduces the first nb rows and columns of a symmetric/Hermitian matrix A to real tridiago
nal form by an unitary similarity transformation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLATRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   E( * ) */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLATRD reduces NB rows and columns of a complex Hermitian matrix A to */
/* > Hermitian tridiagonal form by a unitary similarity */
/* > transformation Q**H * A * Q, and returns the matrices V and W which are */
/* > needed to apply the transformation to the unreduced part of A. */
/* > */
/* > If UPLO = 'U', ZLATRD reduces the last NB rows and columns of a */
/* > matrix, of which the upper triangle is supplied; */
/* > if UPLO = 'L', ZLATRD reduces the first NB rows and columns of a */
/* > matrix, of which the lower triangle is supplied. */
/* > */
/* > This is an auxiliary routine called by ZHETRD. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          Hermitian matrix A is stored: */
/* >          = 'U': Upper triangular */
/* >          = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The number of rows and columns to be reduced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
/* >          n-by-n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n-by-n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit: */
/* >          if UPLO = 'U', the last NB columns have been reduced to */
/* >            tridiagonal form, with the diagonal elements overwriting */
/* >            the diagonal elements of A; the elements above the diagonal */
/* >            with the array TAU, represent the unitary matrix Q as a */
/* >            product of elementary reflectors; */
/* >          if UPLO = 'L', the first NB columns have been reduced to */
/* >            tridiagonal form, with the diagonal elements overwriting */
/* >            the diagonal elements of A; the elements below the diagonal */
/* >            with the array TAU, represent the  unitary matrix Q as a */
/* >            product of elementary reflectors. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal */
/* >          elements of the last NB columns of the reduced matrix; */
/* >          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of */
/* >          the first NB columns of the reduced matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors, stored in */
/* >          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX*16 array, dimension (LDW,NB) */
/* >          The n-by-nb matrix W required to update the unreduced part */
/* >          of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDW */
/* > \verbatim */
/* >          LDW is INTEGER */
/* >          The leading dimension of the array W. LDW >= max(1,N). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(n) H(n-1) . . . H(n-nb+1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i), */
/* >  and tau in TAU(i-1). */
/* > */
/* >  If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(nb). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i), */
/* >  and tau in TAU(i). */
/* > */
/* >  The elements of the vectors v together form the n-by-nb matrix V */
/* >  which is needed, with W, to apply the transformation to the unreduced */
/* >  part of the matrix, using a Hermitian rank-2k update of the form: */
/* >  A := A - V*W**H - W*V**H. */
/* > */
/* >  The contents of A on exit are illustrated by the following examples */
/* >  with n = 5 and nb = 2: */
/* > */
/* >  if UPLO = 'U':                       if UPLO = 'L': */
/* > */
/* >    (  a   a   a   v4  v5 )              (  d                  ) */
/* >    (      a   a   v4  v5 )              (  1   d              ) */
/* >    (          a   1   v5 )              (  v1  1   a          ) */
/* >    (              d   1  )              (  v1  v2  a   a      ) */
/* >    (                  d  )              (  v1  v2  a   a   a  ) */
/* > */
/* >  where d denotes a diagonal element of the reduced matrix, a denotes */
/* >  an element of the original matrix that is unchanged, and vi denotes */
/* >  an element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlatrd_(char *uplo, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *e, doublecomplex *tau, 
	doublecomplex *w, integer *ldw, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Local variables */
    static integer i__, iw;
    static doublecomplex alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zhemv_(char *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), zaxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), zlarfg_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), zlacgv_(integer *, doublecomplex *, 
	    integer *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 243 "zlatrd.f"
    /* Parameter adjustments */
#line 243 "zlatrd.f"
    a_dim1 = *lda;
#line 243 "zlatrd.f"
    a_offset = 1 + a_dim1;
#line 243 "zlatrd.f"
    a -= a_offset;
#line 243 "zlatrd.f"
    --e;
#line 243 "zlatrd.f"
    --tau;
#line 243 "zlatrd.f"
    w_dim1 = *ldw;
#line 243 "zlatrd.f"
    w_offset = 1 + w_dim1;
#line 243 "zlatrd.f"
    w -= w_offset;
#line 243 "zlatrd.f"

#line 243 "zlatrd.f"
    /* Function Body */
#line 243 "zlatrd.f"
    if (*n <= 0) {
#line 243 "zlatrd.f"
	return 0;
#line 243 "zlatrd.f"
    }

#line 246 "zlatrd.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Reduce last NB columns of upper triangle */

#line 250 "zlatrd.f"
	i__1 = *n - *nb + 1;
#line 250 "zlatrd.f"
	for (i__ = *n; i__ >= i__1; --i__) {
#line 251 "zlatrd.f"
	    iw = i__ - *n + *nb;
#line 252 "zlatrd.f"
	    if (i__ < *n) {

/*              Update A(1:i,i) */

#line 256 "zlatrd.f"
		i__2 = i__ + i__ * a_dim1;
#line 256 "zlatrd.f"
		i__3 = i__ + i__ * a_dim1;
#line 256 "zlatrd.f"
		d__1 = a[i__3].r;
#line 256 "zlatrd.f"
		a[i__2].r = d__1, a[i__2].i = 0.;
#line 257 "zlatrd.f"
		i__2 = *n - i__;
#line 257 "zlatrd.f"
		zlacgv_(&i__2, &w[i__ + (iw + 1) * w_dim1], ldw);
#line 258 "zlatrd.f"
		i__2 = *n - i__;
#line 258 "zlatrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 258 "zlatrd.f"
		zgemv_("No transpose", &i__, &i__2, &z__1, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &w[i__ + (iw + 1) * w_dim1], ldw, &
			c_b2, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)12);
#line 260 "zlatrd.f"
		i__2 = *n - i__;
#line 260 "zlatrd.f"
		zlacgv_(&i__2, &w[i__ + (iw + 1) * w_dim1], ldw);
#line 261 "zlatrd.f"
		i__2 = *n - i__;
#line 261 "zlatrd.f"
		zlacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 262 "zlatrd.f"
		i__2 = *n - i__;
#line 262 "zlatrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 262 "zlatrd.f"
		zgemv_("No transpose", &i__, &i__2, &z__1, &w[(iw + 1) * 
			w_dim1 + 1], ldw, &a[i__ + (i__ + 1) * a_dim1], lda, &
			c_b2, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)12);
#line 264 "zlatrd.f"
		i__2 = *n - i__;
#line 264 "zlatrd.f"
		zlacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 265 "zlatrd.f"
		i__2 = i__ + i__ * a_dim1;
#line 265 "zlatrd.f"
		i__3 = i__ + i__ * a_dim1;
#line 265 "zlatrd.f"
		d__1 = a[i__3].r;
#line 265 "zlatrd.f"
		a[i__2].r = d__1, a[i__2].i = 0.;
#line 266 "zlatrd.f"
	    }
#line 267 "zlatrd.f"
	    if (i__ > 1) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(1:i-2,i) */

#line 272 "zlatrd.f"
		i__2 = i__ - 1 + i__ * a_dim1;
#line 272 "zlatrd.f"
		alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 273 "zlatrd.f"
		i__2 = i__ - 1;
#line 273 "zlatrd.f"
		zlarfg_(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, &tau[i__ 
			- 1]);
#line 274 "zlatrd.f"
		i__2 = i__ - 1;
#line 274 "zlatrd.f"
		e[i__2] = alpha.r;
#line 275 "zlatrd.f"
		i__2 = i__ - 1 + i__ * a_dim1;
#line 275 "zlatrd.f"
		a[i__2].r = 1., a[i__2].i = 0.;

/*              Compute W(1:i-1,i) */

#line 279 "zlatrd.f"
		i__2 = i__ - 1;
#line 279 "zlatrd.f"
		zhemv_("Upper", &i__2, &c_b2, &a[a_offset], lda, &a[i__ * 
			a_dim1 + 1], &c__1, &c_b1, &w[iw * w_dim1 + 1], &c__1,
			 (ftnlen)5);
#line 281 "zlatrd.f"
		if (i__ < *n) {
#line 282 "zlatrd.f"
		    i__2 = i__ - 1;
#line 282 "zlatrd.f"
		    i__3 = *n - i__;
#line 282 "zlatrd.f"
		    zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &w[(iw 
			    + 1) * w_dim1 + 1], ldw, &a[i__ * a_dim1 + 1], &
			    c__1, &c_b1, &w[i__ + 1 + iw * w_dim1], &c__1, (
			    ftnlen)19);
#line 285 "zlatrd.f"
		    i__2 = i__ - 1;
#line 285 "zlatrd.f"
		    i__3 = *n - i__;
#line 285 "zlatrd.f"
		    z__1.r = -1., z__1.i = -0.;
#line 285 "zlatrd.f"
		    zgemv_("No transpose", &i__2, &i__3, &z__1, &a[(i__ + 1) *
			     a_dim1 + 1], lda, &w[i__ + 1 + iw * w_dim1], &
			    c__1, &c_b2, &w[iw * w_dim1 + 1], &c__1, (ftnlen)
			    12);
#line 288 "zlatrd.f"
		    i__2 = i__ - 1;
#line 288 "zlatrd.f"
		    i__3 = *n - i__;
#line 288 "zlatrd.f"
		    zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[(
			    i__ + 1) * a_dim1 + 1], lda, &a[i__ * a_dim1 + 1],
			     &c__1, &c_b1, &w[i__ + 1 + iw * w_dim1], &c__1, (
			    ftnlen)19);
#line 291 "zlatrd.f"
		    i__2 = i__ - 1;
#line 291 "zlatrd.f"
		    i__3 = *n - i__;
#line 291 "zlatrd.f"
		    z__1.r = -1., z__1.i = -0.;
#line 291 "zlatrd.f"
		    zgemv_("No transpose", &i__2, &i__3, &z__1, &w[(iw + 1) * 
			    w_dim1 + 1], ldw, &w[i__ + 1 + iw * w_dim1], &
			    c__1, &c_b2, &w[iw * w_dim1 + 1], &c__1, (ftnlen)
			    12);
#line 294 "zlatrd.f"
		}
#line 295 "zlatrd.f"
		i__2 = i__ - 1;
#line 295 "zlatrd.f"
		zscal_(&i__2, &tau[i__ - 1], &w[iw * w_dim1 + 1], &c__1);
#line 296 "zlatrd.f"
		z__3.r = -.5, z__3.i = -0.;
#line 296 "zlatrd.f"
		i__2 = i__ - 1;
#line 296 "zlatrd.f"
		z__2.r = z__3.r * tau[i__2].r - z__3.i * tau[i__2].i, z__2.i =
			 z__3.r * tau[i__2].i + z__3.i * tau[i__2].r;
#line 296 "zlatrd.f"
		i__3 = i__ - 1;
#line 296 "zlatrd.f"
		zdotc_(&z__4, &i__3, &w[iw * w_dim1 + 1], &c__1, &a[i__ * 
			a_dim1 + 1], &c__1);
#line 296 "zlatrd.f"
		z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * 
			z__4.i + z__2.i * z__4.r;
#line 296 "zlatrd.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 298 "zlatrd.f"
		i__2 = i__ - 1;
#line 298 "zlatrd.f"
		zaxpy_(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, &w[iw * 
			w_dim1 + 1], &c__1);
#line 299 "zlatrd.f"
	    }

#line 301 "zlatrd.f"
/* L10: */
#line 301 "zlatrd.f"
	}
#line 302 "zlatrd.f"
    } else {

/*        Reduce first NB columns of lower triangle */

#line 306 "zlatrd.f"
	i__1 = *nb;
#line 306 "zlatrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Update A(i:n,i) */

#line 310 "zlatrd.f"
	    i__2 = i__ + i__ * a_dim1;
#line 310 "zlatrd.f"
	    i__3 = i__ + i__ * a_dim1;
#line 310 "zlatrd.f"
	    d__1 = a[i__3].r;
#line 310 "zlatrd.f"
	    a[i__2].r = d__1, a[i__2].i = 0.;
#line 311 "zlatrd.f"
	    i__2 = i__ - 1;
#line 311 "zlatrd.f"
	    zlacgv_(&i__2, &w[i__ + w_dim1], ldw);
#line 312 "zlatrd.f"
	    i__2 = *n - i__ + 1;
#line 312 "zlatrd.f"
	    i__3 = i__ - 1;
#line 312 "zlatrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 312 "zlatrd.f"
	    zgemv_("No transpose", &i__2, &i__3, &z__1, &a[i__ + a_dim1], lda,
		     &w[i__ + w_dim1], ldw, &c_b2, &a[i__ + i__ * a_dim1], &
		    c__1, (ftnlen)12);
#line 314 "zlatrd.f"
	    i__2 = i__ - 1;
#line 314 "zlatrd.f"
	    zlacgv_(&i__2, &w[i__ + w_dim1], ldw);
#line 315 "zlatrd.f"
	    i__2 = i__ - 1;
#line 315 "zlatrd.f"
	    zlacgv_(&i__2, &a[i__ + a_dim1], lda);
#line 316 "zlatrd.f"
	    i__2 = *n - i__ + 1;
#line 316 "zlatrd.f"
	    i__3 = i__ - 1;
#line 316 "zlatrd.f"
	    z__1.r = -1., z__1.i = -0.;
#line 316 "zlatrd.f"
	    zgemv_("No transpose", &i__2, &i__3, &z__1, &w[i__ + w_dim1], ldw,
		     &a[i__ + a_dim1], lda, &c_b2, &a[i__ + i__ * a_dim1], &
		    c__1, (ftnlen)12);
#line 318 "zlatrd.f"
	    i__2 = i__ - 1;
#line 318 "zlatrd.f"
	    zlacgv_(&i__2, &a[i__ + a_dim1], lda);
#line 319 "zlatrd.f"
	    i__2 = i__ + i__ * a_dim1;
#line 319 "zlatrd.f"
	    i__3 = i__ + i__ * a_dim1;
#line 319 "zlatrd.f"
	    d__1 = a[i__3].r;
#line 319 "zlatrd.f"
	    a[i__2].r = d__1, a[i__2].i = 0.;
#line 320 "zlatrd.f"
	    if (i__ < *n) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(i+2:n,i) */

#line 325 "zlatrd.f"
		i__2 = i__ + 1 + i__ * a_dim1;
#line 325 "zlatrd.f"
		alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 326 "zlatrd.f"
		i__2 = *n - i__;
/* Computing MIN */
#line 326 "zlatrd.f"
		i__3 = i__ + 2;
#line 326 "zlatrd.f"
		zlarfg_(&i__2, &alpha, &a[min(i__3,*n) + i__ * a_dim1], &c__1,
			 &tau[i__]);
#line 328 "zlatrd.f"
		i__2 = i__;
#line 328 "zlatrd.f"
		e[i__2] = alpha.r;
#line 329 "zlatrd.f"
		i__2 = i__ + 1 + i__ * a_dim1;
#line 329 "zlatrd.f"
		a[i__2].r = 1., a[i__2].i = 0.;

/*              Compute W(i+1:n,i) */

#line 333 "zlatrd.f"
		i__2 = *n - i__;
#line 333 "zlatrd.f"
		zhemv_("Lower", &i__2, &c_b2, &a[i__ + 1 + (i__ + 1) * a_dim1]
			, lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b1, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)5);
#line 335 "zlatrd.f"
		i__2 = *n - i__;
#line 335 "zlatrd.f"
		i__3 = i__ - 1;
#line 335 "zlatrd.f"
		zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &w[i__ + 1 
			+ w_dim1], ldw, &a[i__ + 1 + i__ * a_dim1], &c__1, &
			c_b1, &w[i__ * w_dim1 + 1], &c__1, (ftnlen)19);
#line 338 "zlatrd.f"
		i__2 = *n - i__;
#line 338 "zlatrd.f"
		i__3 = i__ - 1;
#line 338 "zlatrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 338 "zlatrd.f"
		zgemv_("No transpose", &i__2, &i__3, &z__1, &a[i__ + 1 + 
			a_dim1], lda, &w[i__ * w_dim1 + 1], &c__1, &c_b2, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)12);
#line 340 "zlatrd.f"
		i__2 = *n - i__;
#line 340 "zlatrd.f"
		i__3 = i__ - 1;
#line 340 "zlatrd.f"
		zgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &a[i__ + 1 
			+ a_dim1], lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &
			c_b1, &w[i__ * w_dim1 + 1], &c__1, (ftnlen)19);
#line 343 "zlatrd.f"
		i__2 = *n - i__;
#line 343 "zlatrd.f"
		i__3 = i__ - 1;
#line 343 "zlatrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 343 "zlatrd.f"
		zgemv_("No transpose", &i__2, &i__3, &z__1, &w[i__ + 1 + 
			w_dim1], ldw, &w[i__ * w_dim1 + 1], &c__1, &c_b2, &w[
			i__ + 1 + i__ * w_dim1], &c__1, (ftnlen)12);
#line 345 "zlatrd.f"
		i__2 = *n - i__;
#line 345 "zlatrd.f"
		zscal_(&i__2, &tau[i__], &w[i__ + 1 + i__ * w_dim1], &c__1);
#line 346 "zlatrd.f"
		z__3.r = -.5, z__3.i = -0.;
#line 346 "zlatrd.f"
		i__2 = i__;
#line 346 "zlatrd.f"
		z__2.r = z__3.r * tau[i__2].r - z__3.i * tau[i__2].i, z__2.i =
			 z__3.r * tau[i__2].i + z__3.i * tau[i__2].r;
#line 346 "zlatrd.f"
		i__3 = *n - i__;
#line 346 "zlatrd.f"
		zdotc_(&z__4, &i__3, &w[i__ + 1 + i__ * w_dim1], &c__1, &a[
			i__ + 1 + i__ * a_dim1], &c__1);
#line 346 "zlatrd.f"
		z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * 
			z__4.i + z__2.i * z__4.r;
#line 346 "zlatrd.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 348 "zlatrd.f"
		i__2 = *n - i__;
#line 348 "zlatrd.f"
		zaxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &w[
			i__ + 1 + i__ * w_dim1], &c__1);
#line 349 "zlatrd.f"
	    }

#line 351 "zlatrd.f"
/* L20: */
#line 351 "zlatrd.f"
	}
#line 352 "zlatrd.f"
    }

#line 354 "zlatrd.f"
    return 0;

/*     End of ZLATRD */

} /* zlatrd_ */

