#line 1 "dsptrd.f"
/* dsptrd.f -- translated by f2c (version 20100827).
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

#line 1 "dsptrd.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = 0.;
static doublereal c_b14 = -1.;

/* > \brief \b DSPTRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPTRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsptrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsptrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsptrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPTRD( UPLO, N, AP, D, E, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), D( * ), E( * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPTRD reduces a real symmetric matrix A stored in packed form to */
/* > symmetric tridiagonal form T by an orthogonal similarity */
/* > transformation: Q**T * A * Q = T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

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
/* >  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP, */
/* >  overwriting A(1:i-1,i+1), and tau is stored in TAU(i). */
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
/* >  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP, */
/* >  overwriting A(i+2:n,i), and tau is stored in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsptrd_(char *uplo, integer *n, doublereal *ap, 
	doublereal *d__, doublereal *e, doublereal *tau, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, i1, ii, i1i1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal taui;
    extern /* Subroutine */ int dspr2_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dspmv_(char *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 190 "dsptrd.f"
    /* Parameter adjustments */
#line 190 "dsptrd.f"
    --tau;
#line 190 "dsptrd.f"
    --e;
#line 190 "dsptrd.f"
    --d__;
#line 190 "dsptrd.f"
    --ap;
#line 190 "dsptrd.f"

#line 190 "dsptrd.f"
    /* Function Body */
#line 190 "dsptrd.f"
    *info = 0;
#line 191 "dsptrd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 192 "dsptrd.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 193 "dsptrd.f"
	*info = -1;
#line 194 "dsptrd.f"
    } else if (*n < 0) {
#line 195 "dsptrd.f"
	*info = -2;
#line 196 "dsptrd.f"
    }
#line 197 "dsptrd.f"
    if (*info != 0) {
#line 198 "dsptrd.f"
	i__1 = -(*info);
#line 198 "dsptrd.f"
	xerbla_("DSPTRD", &i__1, (ftnlen)6);
#line 199 "dsptrd.f"
	return 0;
#line 200 "dsptrd.f"
    }

/*     Quick return if possible */

#line 204 "dsptrd.f"
    if (*n <= 0) {
#line 204 "dsptrd.f"
	return 0;
#line 204 "dsptrd.f"
    }

#line 207 "dsptrd.f"
    if (upper) {

/*        Reduce the upper triangle of A. */
/*        I1 is the index in AP of A(1,I+1). */

#line 212 "dsptrd.f"
	i1 = *n * (*n - 1) / 2 + 1;
#line 213 "dsptrd.f"
	for (i__ = *n - 1; i__ >= 1; --i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**T */
/*           to annihilate A(1:i-1,i+1) */

#line 218 "dsptrd.f"
	    dlarfg_(&i__, &ap[i1 + i__ - 1], &ap[i1], &c__1, &taui);
#line 219 "dsptrd.f"
	    e[i__] = ap[i1 + i__ - 1];

#line 221 "dsptrd.f"
	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

#line 225 "dsptrd.f"
		ap[i1 + i__ - 1] = 1.;

/*              Compute  y := tau * A * v  storing y in TAU(1:i) */

#line 229 "dsptrd.f"
		dspmv_(uplo, &i__, &taui, &ap[1], &ap[i1], &c__1, &c_b8, &tau[
			1], &c__1, (ftnlen)1);

/*              Compute  w := y - 1/2 * tau * (y**T *v) * v */

#line 234 "dsptrd.f"
		alpha = taui * -.5 * ddot_(&i__, &tau[1], &c__1, &ap[i1], &
			c__1);
#line 235 "dsptrd.f"
		daxpy_(&i__, &alpha, &ap[i1], &c__1, &tau[1], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**T - w * v**T */

#line 240 "dsptrd.f"
		dspr2_(uplo, &i__, &c_b14, &ap[i1], &c__1, &tau[1], &c__1, &
			ap[1], (ftnlen)1);

#line 242 "dsptrd.f"
		ap[i1 + i__ - 1] = e[i__];
#line 243 "dsptrd.f"
	    }
#line 244 "dsptrd.f"
	    d__[i__ + 1] = ap[i1 + i__];
#line 245 "dsptrd.f"
	    tau[i__] = taui;
#line 246 "dsptrd.f"
	    i1 -= i__;
#line 247 "dsptrd.f"
/* L10: */
#line 247 "dsptrd.f"
	}
#line 248 "dsptrd.f"
	d__[1] = ap[1];
#line 249 "dsptrd.f"
    } else {

/*        Reduce the lower triangle of A. II is the index in AP of */
/*        A(i,i) and I1I1 is the index of A(i+1,i+1). */

#line 254 "dsptrd.f"
	ii = 1;
#line 255 "dsptrd.f"
	i__1 = *n - 1;
#line 255 "dsptrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 256 "dsptrd.f"
	    i1i1 = ii + *n - i__ + 1;

/*           Generate elementary reflector H(i) = I - tau * v * v**T */
/*           to annihilate A(i+2:n,i) */

#line 261 "dsptrd.f"
	    i__2 = *n - i__;
#line 261 "dsptrd.f"
	    dlarfg_(&i__2, &ap[ii + 1], &ap[ii + 2], &c__1, &taui);
#line 262 "dsptrd.f"
	    e[i__] = ap[ii + 1];

#line 264 "dsptrd.f"
	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) */

#line 268 "dsptrd.f"
		ap[ii + 1] = 1.;

/*              Compute  y := tau * A * v  storing y in TAU(i:n-1) */

#line 272 "dsptrd.f"
		i__2 = *n - i__;
#line 272 "dsptrd.f"
		dspmv_(uplo, &i__2, &taui, &ap[i1i1], &ap[ii + 1], &c__1, &
			c_b8, &tau[i__], &c__1, (ftnlen)1);

/*              Compute  w := y - 1/2 * tau * (y**T *v) * v */

#line 277 "dsptrd.f"
		i__2 = *n - i__;
#line 277 "dsptrd.f"
		alpha = taui * -.5 * ddot_(&i__2, &tau[i__], &c__1, &ap[ii + 
			1], &c__1);
#line 279 "dsptrd.f"
		i__2 = *n - i__;
#line 279 "dsptrd.f"
		daxpy_(&i__2, &alpha, &ap[ii + 1], &c__1, &tau[i__], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**T - w * v**T */

#line 284 "dsptrd.f"
		i__2 = *n - i__;
#line 284 "dsptrd.f"
		dspr2_(uplo, &i__2, &c_b14, &ap[ii + 1], &c__1, &tau[i__], &
			c__1, &ap[i1i1], (ftnlen)1);

#line 287 "dsptrd.f"
		ap[ii + 1] = e[i__];
#line 288 "dsptrd.f"
	    }
#line 289 "dsptrd.f"
	    d__[i__] = ap[ii];
#line 290 "dsptrd.f"
	    tau[i__] = taui;
#line 291 "dsptrd.f"
	    ii = i1i1;
#line 292 "dsptrd.f"
/* L20: */
#line 292 "dsptrd.f"
	}
#line 293 "dsptrd.f"
	d__[*n] = ap[ii];
#line 294 "dsptrd.f"
    }

#line 296 "dsptrd.f"
    return 0;

/*     End of DSPTRD */

} /* dsptrd_ */

