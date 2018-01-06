#line 1 "chptrd.f"
/* chptrd.f -- translated by f2c (version 20100827).
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

#line 1 "chptrd.f"
/* Table of constant values */

static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b CHPTRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPTRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chptrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chptrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chptrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPTRD( UPLO, N, AP, D, E, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ) */
/*       COMPLEX            AP( * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPTRD reduces a complex Hermitian matrix A stored in packed form to */
/* > real symmetric tridiagonal form T by a unitary similarity */
/* > transformation: Q**H * A * Q = T. */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* >          On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* >          of A are overwritten by the corresponding elements of the */
/* >          tridiagonal matrix T, and the elements above the first */
/* >          superdiagonal, with the array TAU, represent the unitary */
/* >          matrix Q as a product of elementary reflectors; if UPLO */
/* >          = 'L', the diagonal and first subdiagonal of A are over- */
/* >          written by the corresponding elements of the tridiagonal */
/* >          matrix T, and the elements below the first subdiagonal, with */
/* >          the array TAU, represent the unitary matrix Q as a product */
/* >          of elementary reflectors. See Further Details. */
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
/* >          TAU is COMPLEX array, dimension (N-1) */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

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
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
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
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP, */
/* >  overwriting A(i+2:n,i), and tau is stored in TAU(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int chptrd_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *d__, doublereal *e, doublecomplex *tau, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Local variables */
    static integer i__, i1, ii, i1i1;
    static doublecomplex taui;
    extern /* Subroutine */ int chpr2_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, ftnlen);
    static doublecomplex alpha;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int chpmv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), caxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int clarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *), xerbla_(char *, 
	    integer *, ftnlen);


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

#line 196 "chptrd.f"
    /* Parameter adjustments */
#line 196 "chptrd.f"
    --tau;
#line 196 "chptrd.f"
    --e;
#line 196 "chptrd.f"
    --d__;
#line 196 "chptrd.f"
    --ap;
#line 196 "chptrd.f"

#line 196 "chptrd.f"
    /* Function Body */
#line 196 "chptrd.f"
    *info = 0;
#line 197 "chptrd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 198 "chptrd.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 199 "chptrd.f"
	*info = -1;
#line 200 "chptrd.f"
    } else if (*n < 0) {
#line 201 "chptrd.f"
	*info = -2;
#line 202 "chptrd.f"
    }
#line 203 "chptrd.f"
    if (*info != 0) {
#line 204 "chptrd.f"
	i__1 = -(*info);
#line 204 "chptrd.f"
	xerbla_("CHPTRD", &i__1, (ftnlen)6);
#line 205 "chptrd.f"
	return 0;
#line 206 "chptrd.f"
    }

/*     Quick return if possible */

#line 210 "chptrd.f"
    if (*n <= 0) {
#line 210 "chptrd.f"
	return 0;
#line 210 "chptrd.f"
    }

#line 213 "chptrd.f"
    if (upper) {

/*        Reduce the upper triangle of A. */
/*        I1 is the index in AP of A(1,I+1). */

#line 218 "chptrd.f"
	i1 = *n * (*n - 1) / 2 + 1;
#line 219 "chptrd.f"
	i__1 = i1 + *n - 1;
#line 219 "chptrd.f"
	i__2 = i1 + *n - 1;
#line 219 "chptrd.f"
	d__1 = ap[i__2].r;
#line 219 "chptrd.f"
	ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 220 "chptrd.f"
	for (i__ = *n - 1; i__ >= 1; --i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**H */
/*           to annihilate A(1:i-1,i+1) */

#line 225 "chptrd.f"
	    i__1 = i1 + i__ - 1;
#line 225 "chptrd.f"
	    alpha.r = ap[i__1].r, alpha.i = ap[i__1].i;
#line 226 "chptrd.f"
	    clarfg_(&i__, &alpha, &ap[i1], &c__1, &taui);
#line 227 "chptrd.f"
	    i__1 = i__;
#line 227 "chptrd.f"
	    e[i__1] = alpha.r;

#line 229 "chptrd.f"
	    if (taui.r != 0. || taui.i != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

#line 233 "chptrd.f"
		i__1 = i1 + i__ - 1;
#line 233 "chptrd.f"
		ap[i__1].r = 1., ap[i__1].i = 0.;

/*              Compute  y := tau * A * v  storing y in TAU(1:i) */

#line 237 "chptrd.f"
		chpmv_(uplo, &i__, &taui, &ap[1], &ap[i1], &c__1, &c_b2, &tau[
			1], &c__1, (ftnlen)1);

/*              Compute  w := y - 1/2 * tau * (y**H *v) * v */

#line 242 "chptrd.f"
		z__3.r = -.5, z__3.i = -0.;
#line 242 "chptrd.f"
		z__2.r = z__3.r * taui.r - z__3.i * taui.i, z__2.i = z__3.r * 
			taui.i + z__3.i * taui.r;
#line 242 "chptrd.f"
		cdotc_(&z__4, &i__, &tau[1], &c__1, &ap[i1], &c__1);
#line 242 "chptrd.f"
		z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * 
			z__4.i + z__2.i * z__4.r;
#line 242 "chptrd.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 243 "chptrd.f"
		caxpy_(&i__, &alpha, &ap[i1], &c__1, &tau[1], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**H - w * v**H */

#line 248 "chptrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 248 "chptrd.f"
		chpr2_(uplo, &i__, &z__1, &ap[i1], &c__1, &tau[1], &c__1, &ap[
			1], (ftnlen)1);

#line 250 "chptrd.f"
	    }
#line 251 "chptrd.f"
	    i__1 = i1 + i__ - 1;
#line 251 "chptrd.f"
	    i__2 = i__;
#line 251 "chptrd.f"
	    ap[i__1].r = e[i__2], ap[i__1].i = 0.;
#line 252 "chptrd.f"
	    i__1 = i__ + 1;
#line 252 "chptrd.f"
	    i__2 = i1 + i__;
#line 252 "chptrd.f"
	    d__[i__1] = ap[i__2].r;
#line 253 "chptrd.f"
	    i__1 = i__;
#line 253 "chptrd.f"
	    tau[i__1].r = taui.r, tau[i__1].i = taui.i;
#line 254 "chptrd.f"
	    i1 -= i__;
#line 255 "chptrd.f"
/* L10: */
#line 255 "chptrd.f"
	}
#line 256 "chptrd.f"
	d__[1] = ap[1].r;
#line 257 "chptrd.f"
    } else {

/*        Reduce the lower triangle of A. II is the index in AP of */
/*        A(i,i) and I1I1 is the index of A(i+1,i+1). */

#line 262 "chptrd.f"
	ii = 1;
#line 263 "chptrd.f"
	d__1 = ap[1].r;
#line 263 "chptrd.f"
	ap[1].r = d__1, ap[1].i = 0.;
#line 264 "chptrd.f"
	i__1 = *n - 1;
#line 264 "chptrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 265 "chptrd.f"
	    i1i1 = ii + *n - i__ + 1;

/*           Generate elementary reflector H(i) = I - tau * v * v**H */
/*           to annihilate A(i+2:n,i) */

#line 270 "chptrd.f"
	    i__2 = ii + 1;
#line 270 "chptrd.f"
	    alpha.r = ap[i__2].r, alpha.i = ap[i__2].i;
#line 271 "chptrd.f"
	    i__2 = *n - i__;
#line 271 "chptrd.f"
	    clarfg_(&i__2, &alpha, &ap[ii + 2], &c__1, &taui);
#line 272 "chptrd.f"
	    i__2 = i__;
#line 272 "chptrd.f"
	    e[i__2] = alpha.r;

#line 274 "chptrd.f"
	    if (taui.r != 0. || taui.i != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) */

#line 278 "chptrd.f"
		i__2 = ii + 1;
#line 278 "chptrd.f"
		ap[i__2].r = 1., ap[i__2].i = 0.;

/*              Compute  y := tau * A * v  storing y in TAU(i:n-1) */

#line 282 "chptrd.f"
		i__2 = *n - i__;
#line 282 "chptrd.f"
		chpmv_(uplo, &i__2, &taui, &ap[i1i1], &ap[ii + 1], &c__1, &
			c_b2, &tau[i__], &c__1, (ftnlen)1);

/*              Compute  w := y - 1/2 * tau * (y**H *v) * v */

#line 287 "chptrd.f"
		z__3.r = -.5, z__3.i = -0.;
#line 287 "chptrd.f"
		z__2.r = z__3.r * taui.r - z__3.i * taui.i, z__2.i = z__3.r * 
			taui.i + z__3.i * taui.r;
#line 287 "chptrd.f"
		i__2 = *n - i__;
#line 287 "chptrd.f"
		cdotc_(&z__4, &i__2, &tau[i__], &c__1, &ap[ii + 1], &c__1);
#line 287 "chptrd.f"
		z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r * 
			z__4.i + z__2.i * z__4.r;
#line 287 "chptrd.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 289 "chptrd.f"
		i__2 = *n - i__;
#line 289 "chptrd.f"
		caxpy_(&i__2, &alpha, &ap[ii + 1], &c__1, &tau[i__], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**H - w * v**H */

#line 294 "chptrd.f"
		i__2 = *n - i__;
#line 294 "chptrd.f"
		z__1.r = -1., z__1.i = -0.;
#line 294 "chptrd.f"
		chpr2_(uplo, &i__2, &z__1, &ap[ii + 1], &c__1, &tau[i__], &
			c__1, &ap[i1i1], (ftnlen)1);

#line 297 "chptrd.f"
	    }
#line 298 "chptrd.f"
	    i__2 = ii + 1;
#line 298 "chptrd.f"
	    i__3 = i__;
#line 298 "chptrd.f"
	    ap[i__2].r = e[i__3], ap[i__2].i = 0.;
#line 299 "chptrd.f"
	    i__2 = i__;
#line 299 "chptrd.f"
	    i__3 = ii;
#line 299 "chptrd.f"
	    d__[i__2] = ap[i__3].r;
#line 300 "chptrd.f"
	    i__2 = i__;
#line 300 "chptrd.f"
	    tau[i__2].r = taui.r, tau[i__2].i = taui.i;
#line 301 "chptrd.f"
	    ii = i1i1;
#line 302 "chptrd.f"
/* L20: */
#line 302 "chptrd.f"
	}
#line 303 "chptrd.f"
	i__1 = *n;
#line 303 "chptrd.f"
	i__2 = ii;
#line 303 "chptrd.f"
	d__[i__1] = ap[i__2].r;
#line 304 "chptrd.f"
    }

#line 306 "chptrd.f"
    return 0;

/*     End of CHPTRD */

} /* chptrd_ */

