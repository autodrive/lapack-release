#line 1 "ssytd2.f"
/* ssytd2.f -- translated by f2c (version 20100827).
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

#line 1 "ssytd2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = 0.;
static doublereal c_b14 = -1.;

/* > \brief \b SSYTD2 reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarit
y transformation (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytd2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytd2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytd2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), D( * ), E( * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          D is REAL array, dimension (N) */
/* >          The diagonal elements of the tridiagonal matrix T: */
/* >          D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          The off-diagonal elements of the tridiagonal matrix T: */
/* >          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (N-1) */
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

/* > \ingroup realSYcomputational */

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
/* Subroutine */ int ssytd2_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__;
    static doublereal taui;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int ssyr2_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), ssymv_(char *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), slarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);


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

#line 215 "ssytd2.f"
    /* Parameter adjustments */
#line 215 "ssytd2.f"
    a_dim1 = *lda;
#line 215 "ssytd2.f"
    a_offset = 1 + a_dim1;
#line 215 "ssytd2.f"
    a -= a_offset;
#line 215 "ssytd2.f"
    --d__;
#line 215 "ssytd2.f"
    --e;
#line 215 "ssytd2.f"
    --tau;
#line 215 "ssytd2.f"

#line 215 "ssytd2.f"
    /* Function Body */
#line 215 "ssytd2.f"
    *info = 0;
#line 216 "ssytd2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 217 "ssytd2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 218 "ssytd2.f"
	*info = -1;
#line 219 "ssytd2.f"
    } else if (*n < 0) {
#line 220 "ssytd2.f"
	*info = -2;
#line 221 "ssytd2.f"
    } else if (*lda < max(1,*n)) {
#line 222 "ssytd2.f"
	*info = -4;
#line 223 "ssytd2.f"
    }
#line 224 "ssytd2.f"
    if (*info != 0) {
#line 225 "ssytd2.f"
	i__1 = -(*info);
#line 225 "ssytd2.f"
	xerbla_("SSYTD2", &i__1, (ftnlen)6);
#line 226 "ssytd2.f"
	return 0;
#line 227 "ssytd2.f"
    }

/*     Quick return if possible */

#line 231 "ssytd2.f"
    if (*n <= 0) {
#line 231 "ssytd2.f"
	return 0;
#line 231 "ssytd2.f"
    }

#line 234 "ssytd2.f"
    if (upper) {

/*        Reduce the upper triangle of A */

#line 238 "ssytd2.f"
	for (i__ = *n - 1; i__ >= 1; --i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**T */
/*           to annihilate A(1:i-1,i+1) */

#line 243 "ssytd2.f"
	    slarfg_(&i__, &a[i__ + (i__ + 1) * a_dim1], &a[(i__ + 1) * a_dim1 
		    + 1], &c__1, &taui);
#line 244 "ssytd2.f"
	    e[i__] = a[i__ + (i__ + 1) * a_dim1];

#line 246 "ssytd2.f"
	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

#line 250 "ssytd2.f"
		a[i__ + (i__ + 1) * a_dim1] = 1.;

/*              Compute  x := tau * A * v  storing x in TAU(1:i) */

#line 254 "ssytd2.f"
		ssymv_(uplo, &i__, &taui, &a[a_offset], lda, &a[(i__ + 1) * 
			a_dim1 + 1], &c__1, &c_b8, &tau[1], &c__1, (ftnlen)1);

/*              Compute  w := x - 1/2 * tau * (x**T * v) * v */

#line 259 "ssytd2.f"
		alpha = taui * -.5 * sdot_(&i__, &tau[1], &c__1, &a[(i__ + 1) 
			* a_dim1 + 1], &c__1);
#line 260 "ssytd2.f"
		saxpy_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[
			1], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**T - w * v**T */

#line 265 "ssytd2.f"
		ssyr2_(uplo, &i__, &c_b14, &a[(i__ + 1) * a_dim1 + 1], &c__1, 
			&tau[1], &c__1, &a[a_offset], lda, (ftnlen)1);

#line 268 "ssytd2.f"
		a[i__ + (i__ + 1) * a_dim1] = e[i__];
#line 269 "ssytd2.f"
	    }
#line 270 "ssytd2.f"
	    d__[i__ + 1] = a[i__ + 1 + (i__ + 1) * a_dim1];
#line 271 "ssytd2.f"
	    tau[i__] = taui;
#line 272 "ssytd2.f"
/* L10: */
#line 272 "ssytd2.f"
	}
#line 273 "ssytd2.f"
	d__[1] = a[a_dim1 + 1];
#line 274 "ssytd2.f"
    } else {

/*        Reduce the lower triangle of A */

#line 278 "ssytd2.f"
	i__1 = *n - 1;
#line 278 "ssytd2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**T */
/*           to annihilate A(i+2:n,i) */

#line 283 "ssytd2.f"
	    i__2 = *n - i__;
/* Computing MIN */
#line 283 "ssytd2.f"
	    i__3 = i__ + 2;
#line 283 "ssytd2.f"
	    slarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + i__ *
		     a_dim1], &c__1, &taui);
#line 285 "ssytd2.f"
	    e[i__] = a[i__ + 1 + i__ * a_dim1];

#line 287 "ssytd2.f"
	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) */

#line 291 "ssytd2.f"
		a[i__ + 1 + i__ * a_dim1] = 1.;

/*              Compute  x := tau * A * v  storing y in TAU(i:n-1) */

#line 295 "ssytd2.f"
		i__2 = *n - i__;
#line 295 "ssytd2.f"
		ssymv_(uplo, &i__2, &taui, &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b8, &tau[
			i__], &c__1, (ftnlen)1);

/*              Compute  w := x - 1/2 * tau * (x**T * v) * v */

#line 300 "ssytd2.f"
		i__2 = *n - i__;
#line 300 "ssytd2.f"
		alpha = taui * -.5 * sdot_(&i__2, &tau[i__], &c__1, &a[i__ + 
			1 + i__ * a_dim1], &c__1);
#line 302 "ssytd2.f"
		i__2 = *n - i__;
#line 302 "ssytd2.f"
		saxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[
			i__], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**T - w * v**T */

#line 307 "ssytd2.f"
		i__2 = *n - i__;
#line 307 "ssytd2.f"
		ssyr2_(uplo, &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1,
			 &tau[i__], &c__1, &a[i__ + 1 + (i__ + 1) * a_dim1], 
			lda, (ftnlen)1);

#line 310 "ssytd2.f"
		a[i__ + 1 + i__ * a_dim1] = e[i__];
#line 311 "ssytd2.f"
	    }
#line 312 "ssytd2.f"
	    d__[i__] = a[i__ + i__ * a_dim1];
#line 313 "ssytd2.f"
	    tau[i__] = taui;
#line 314 "ssytd2.f"
/* L20: */
#line 314 "ssytd2.f"
	}
#line 315 "ssytd2.f"
	d__[*n] = a[*n + *n * a_dim1];
#line 316 "ssytd2.f"
    }

#line 318 "ssytd2.f"
    return 0;

/*     End of SSYTD2 */

} /* ssytd2_ */

