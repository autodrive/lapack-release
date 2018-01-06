#line 1 "dsytd2.f"
/* dsytd2.f -- translated by f2c (version 20100827).
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

#line 1 "dsytd2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = 0.;
static doublereal c_b14 = -1.;

/* > \brief \b DSYTD2 reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarit
y transformation (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytd2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytd2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytd2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal */
/* > form T by an orthogonal similarity transformation: Q**T * A * Q = T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored: */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          n-by-n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n-by-n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* >          of A are overwritten by the corresponding elements of the */
/* >          tridiagonal matrix T, and the elements above the first */
/* >          superdiagonal, with the array TAU, represent the orthogonal */
/* >          matrix Q as a product of elementary reflectors; if UPLO */
/* >          = 'L', the diagonal and first subdiagonal of A are over- */
/* >          written by the corresponding elements of the tridiagonal */
/* >          matrix T, and the elements below the first subdiagonal, with */
/* >          the array TAU, represent the orthogonal matrix Q as a product */
/* >          of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of the tridiagonal matrix T: */
/* >          D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The off-diagonal elements of the tridiagonal matrix T: */
/* >          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
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

/* > \ingroup doubleSYcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(n-1) . . . H(2) H(1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in */
/* >  A(1:i-1,i+1), and tau in TAU(i). */
/* > */
/* >  If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(n-1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i), */
/* >  and tau in TAU(i). */
/* > */
/* >  The contents of A on exit are illustrated by the following examples */
/* >  with n = 5: */
/* > */
/* >  if UPLO = 'U':                       if UPLO = 'L': */
/* > */
/* >    (  d   e   v2  v3  v4 )              (  d                  ) */
/* >    (      d   e   v3  v4 )              (  e   d              ) */
/* >    (          d   e   v4 )              (  v1  e   d          ) */
/* >    (              d   e  )              (  v1  v2  e   d      ) */
/* >    (                  d  )              (  v1  v2  v3  e   d  ) */
/* > */
/* >  where d and e denote diagonal and off-diagonal elements of T, and vi */
/* >  denotes an element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsytd2_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal taui;
    extern /* Subroutine */ int dsyr2_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlarfg_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *), xerbla_(char *, integer *
	    , ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 216 "dsytd2.f"
    /* Parameter adjustments */
#line 216 "dsytd2.f"
    a_dim1 = *lda;
#line 216 "dsytd2.f"
    a_offset = 1 + a_dim1;
#line 216 "dsytd2.f"
    a -= a_offset;
#line 216 "dsytd2.f"
    --d__;
#line 216 "dsytd2.f"
    --e;
#line 216 "dsytd2.f"
    --tau;
#line 216 "dsytd2.f"

#line 216 "dsytd2.f"
    /* Function Body */
#line 216 "dsytd2.f"
    *info = 0;
#line 217 "dsytd2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 218 "dsytd2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 219 "dsytd2.f"
	*info = -1;
#line 220 "dsytd2.f"
    } else if (*n < 0) {
#line 221 "dsytd2.f"
	*info = -2;
#line 222 "dsytd2.f"
    } else if (*lda < max(1,*n)) {
#line 223 "dsytd2.f"
	*info = -4;
#line 224 "dsytd2.f"
    }
#line 225 "dsytd2.f"
    if (*info != 0) {
#line 226 "dsytd2.f"
	i__1 = -(*info);
#line 226 "dsytd2.f"
	xerbla_("DSYTD2", &i__1, (ftnlen)6);
#line 227 "dsytd2.f"
	return 0;
#line 228 "dsytd2.f"
    }

/*     Quick return if possible */

#line 232 "dsytd2.f"
    if (*n <= 0) {
#line 232 "dsytd2.f"
	return 0;
#line 232 "dsytd2.f"
    }

#line 235 "dsytd2.f"
    if (upper) {

/*        Reduce the upper triangle of A */

#line 239 "dsytd2.f"
	for (i__ = *n - 1; i__ >= 1; --i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**T */
/*           to annihilate A(1:i-1,i+1) */

#line 244 "dsytd2.f"
	    dlarfg_(&i__, &a[i__ + (i__ + 1) * a_dim1], &a[(i__ + 1) * a_dim1 
		    + 1], &c__1, &taui);
#line 245 "dsytd2.f"
	    e[i__] = a[i__ + (i__ + 1) * a_dim1];

#line 247 "dsytd2.f"
	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

#line 251 "dsytd2.f"
		a[i__ + (i__ + 1) * a_dim1] = 1.;

/*              Compute  x := tau * A * v  storing x in TAU(1:i) */

#line 255 "dsytd2.f"
		dsymv_(uplo, &i__, &taui, &a[a_offset], lda, &a[(i__ + 1) * 
			a_dim1 + 1], &c__1, &c_b8, &tau[1], &c__1, (ftnlen)1);

/*              Compute  w := x - 1/2 * tau * (x**T * v) * v */

#line 260 "dsytd2.f"
		alpha = taui * -.5 * ddot_(&i__, &tau[1], &c__1, &a[(i__ + 1) 
			* a_dim1 + 1], &c__1);
#line 261 "dsytd2.f"
		daxpy_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[
			1], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**T - w * v**T */

#line 266 "dsytd2.f"
		dsyr2_(uplo, &i__, &c_b14, &a[(i__ + 1) * a_dim1 + 1], &c__1, 
			&tau[1], &c__1, &a[a_offset], lda, (ftnlen)1);

#line 269 "dsytd2.f"
		a[i__ + (i__ + 1) * a_dim1] = e[i__];
#line 270 "dsytd2.f"
	    }
#line 271 "dsytd2.f"
	    d__[i__ + 1] = a[i__ + 1 + (i__ + 1) * a_dim1];
#line 272 "dsytd2.f"
	    tau[i__] = taui;
#line 273 "dsytd2.f"
/* L10: */
#line 273 "dsytd2.f"
	}
#line 274 "dsytd2.f"
	d__[1] = a[a_dim1 + 1];
#line 275 "dsytd2.f"
    } else {

/*        Reduce the lower triangle of A */

#line 279 "dsytd2.f"
	i__1 = *n - 1;
#line 279 "dsytd2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**T */
/*           to annihilate A(i+2:n,i) */

#line 284 "dsytd2.f"
	    i__2 = *n - i__;
/* Computing MIN */
#line 284 "dsytd2.f"
	    i__3 = i__ + 2;
#line 284 "dsytd2.f"
	    dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + i__ *
		     a_dim1], &c__1, &taui);
#line 286 "dsytd2.f"
	    e[i__] = a[i__ + 1 + i__ * a_dim1];

#line 288 "dsytd2.f"
	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) */

#line 292 "dsytd2.f"
		a[i__ + 1 + i__ * a_dim1] = 1.;

/*              Compute  x := tau * A * v  storing y in TAU(i:n-1) */

#line 296 "dsytd2.f"
		i__2 = *n - i__;
#line 296 "dsytd2.f"
		dsymv_(uplo, &i__2, &taui, &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b8, &tau[
			i__], &c__1, (ftnlen)1);

/*              Compute  w := x - 1/2 * tau * (x**T * v) * v */

#line 301 "dsytd2.f"
		i__2 = *n - i__;
#line 301 "dsytd2.f"
		alpha = taui * -.5 * ddot_(&i__2, &tau[i__], &c__1, &a[i__ + 
			1 + i__ * a_dim1], &c__1);
#line 303 "dsytd2.f"
		i__2 = *n - i__;
#line 303 "dsytd2.f"
		daxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[
			i__], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**T - w * v**T */

#line 308 "dsytd2.f"
		i__2 = *n - i__;
#line 308 "dsytd2.f"
		dsyr2_(uplo, &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1,
			 &tau[i__], &c__1, &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda, (ftnlen)1);

#line 311 "dsytd2.f"
		a[i__ + 1 + i__ * a_dim1] = e[i__];
#line 312 "dsytd2.f"
	    }
#line 313 "dsytd2.f"
	    d__[i__] = a[i__ + i__ * a_dim1];
#line 314 "dsytd2.f"
	    tau[i__] = taui;
#line 315 "dsytd2.f"
/* L20: */
#line 315 "dsytd2.f"
	}
#line 316 "dsytd2.f"
	d__[*n] = a[*n + *n * a_dim1];
#line 317 "dsytd2.f"
    }

#line 319 "dsytd2.f"
    return 0;

/*     End of DSYTD2 */

} /* dsytd2_ */

