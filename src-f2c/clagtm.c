#line 1 "clagtm.f"
/* clagtm.f -- translated by f2c (version 20100827).
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

#line 1 "clagtm.f"
/* > \brief \b CLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matr
ix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAGTM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clagtm.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clagtm.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clagtm.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, */
/*                          B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            LDB, LDX, N, NRHS */
/*       REAL               ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAGTM performs a matrix-vector product of the form */
/* > */
/* >    B := alpha * A * X + beta * B */
/* > */
/* > where A is a tridiagonal matrix of order N, B and X are N by NRHS */
/* > matrices, and alpha and beta are real scalars, each of which may be */
/* > 0., 1., or -1. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the operation applied to A. */
/* >          = 'N':  No transpose, B := alpha * A * X + beta * B */
/* >          = 'T':  Transpose,    B := alpha * A**T * X + beta * B */
/* >          = 'C':  Conjugate transpose, B := alpha * A**H * X + beta * B */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices X and B. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL */
/* >          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise, */
/* >          it is assumed to be 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is COMPLEX array, dimension (N-1) */
/* >          The (n-1) sub-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          The diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is COMPLEX array, dimension (N-1) */
/* >          The (n-1) super-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (LDX,NRHS) */
/* >          The N by NRHS matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X.  LDX >= max(N,1). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is REAL */
/* >          The scalar beta.  BETA must be 0., 1., or -1.; otherwise, */
/* >          it is assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
/* >          On entry, the N by NRHS matrix B. */
/* >          On exit, B is overwritten by the matrix expression */
/* >          B := alpha * A * X + beta * B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(N,1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clagtm_(char *trans, integer *n, integer *nrhs, 
	doublereal *alpha, doublecomplex *dl, doublecomplex *d__, 
	doublecomplex *du, doublecomplex *x, integer *ldx, doublereal *beta, 
	doublecomplex *b, integer *ldb, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9, i__10;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 181 "clagtm.f"
    /* Parameter adjustments */
#line 181 "clagtm.f"
    --dl;
#line 181 "clagtm.f"
    --d__;
#line 181 "clagtm.f"
    --du;
#line 181 "clagtm.f"
    x_dim1 = *ldx;
#line 181 "clagtm.f"
    x_offset = 1 + x_dim1;
#line 181 "clagtm.f"
    x -= x_offset;
#line 181 "clagtm.f"
    b_dim1 = *ldb;
#line 181 "clagtm.f"
    b_offset = 1 + b_dim1;
#line 181 "clagtm.f"
    b -= b_offset;
#line 181 "clagtm.f"

#line 181 "clagtm.f"
    /* Function Body */
#line 181 "clagtm.f"
    if (*n == 0) {
#line 181 "clagtm.f"
	return 0;
#line 181 "clagtm.f"
    }

/*     Multiply B by BETA if BETA.NE.1. */

#line 186 "clagtm.f"
    if (*beta == 0.) {
#line 187 "clagtm.f"
	i__1 = *nrhs;
#line 187 "clagtm.f"
	for (j = 1; j <= i__1; ++j) {
#line 188 "clagtm.f"
	    i__2 = *n;
#line 188 "clagtm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 189 "clagtm.f"
		i__3 = i__ + j * b_dim1;
#line 189 "clagtm.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 190 "clagtm.f"
/* L10: */
#line 190 "clagtm.f"
	    }
#line 191 "clagtm.f"
/* L20: */
#line 191 "clagtm.f"
	}
#line 192 "clagtm.f"
    } else if (*beta == -1.) {
#line 193 "clagtm.f"
	i__1 = *nrhs;
#line 193 "clagtm.f"
	for (j = 1; j <= i__1; ++j) {
#line 194 "clagtm.f"
	    i__2 = *n;
#line 194 "clagtm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 195 "clagtm.f"
		i__3 = i__ + j * b_dim1;
#line 195 "clagtm.f"
		i__4 = i__ + j * b_dim1;
#line 195 "clagtm.f"
		z__1.r = -b[i__4].r, z__1.i = -b[i__4].i;
#line 195 "clagtm.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 196 "clagtm.f"
/* L30: */
#line 196 "clagtm.f"
	    }
#line 197 "clagtm.f"
/* L40: */
#line 197 "clagtm.f"
	}
#line 198 "clagtm.f"
    }

#line 200 "clagtm.f"
    if (*alpha == 1.) {
#line 201 "clagtm.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B + A*X */

#line 205 "clagtm.f"
	    i__1 = *nrhs;
#line 205 "clagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 206 "clagtm.f"
		if (*n == 1) {
#line 207 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 207 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 207 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 207 "clagtm.f"
		    z__2.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i, 
			    z__2.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4]
			    .r;
#line 207 "clagtm.f"
		    z__1.r = b[i__3].r + z__2.r, z__1.i = b[i__3].i + z__2.i;
#line 207 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 208 "clagtm.f"
		} else {
#line 209 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 209 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 209 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 209 "clagtm.f"
		    z__3.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i, 
			    z__3.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4]
			    .r;
#line 209 "clagtm.f"
		    z__2.r = b[i__3].r + z__3.r, z__2.i = b[i__3].i + z__3.i;
#line 209 "clagtm.f"
		    i__5 = j * x_dim1 + 2;
#line 209 "clagtm.f"
		    z__4.r = du[1].r * x[i__5].r - du[1].i * x[i__5].i, 
			    z__4.i = du[1].r * x[i__5].i + du[1].i * x[i__5]
			    .r;
#line 209 "clagtm.f"
		    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 209 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 211 "clagtm.f"
		    i__2 = *n + j * b_dim1;
#line 211 "clagtm.f"
		    i__3 = *n + j * b_dim1;
#line 211 "clagtm.f"
		    i__4 = *n - 1;
#line 211 "clagtm.f"
		    i__5 = *n - 1 + j * x_dim1;
#line 211 "clagtm.f"
		    z__3.r = dl[i__4].r * x[i__5].r - dl[i__4].i * x[i__5].i, 
			    z__3.i = dl[i__4].r * x[i__5].i + dl[i__4].i * x[
			    i__5].r;
#line 211 "clagtm.f"
		    z__2.r = b[i__3].r + z__3.r, z__2.i = b[i__3].i + z__3.i;
#line 211 "clagtm.f"
		    i__6 = *n;
#line 211 "clagtm.f"
		    i__7 = *n + j * x_dim1;
#line 211 "clagtm.f"
		    z__4.r = d__[i__6].r * x[i__7].r - d__[i__6].i * x[i__7]
			    .i, z__4.i = d__[i__6].r * x[i__7].i + d__[i__6]
			    .i * x[i__7].r;
#line 211 "clagtm.f"
		    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 211 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 213 "clagtm.f"
		    i__2 = *n - 1;
#line 213 "clagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 214 "clagtm.f"
			i__3 = i__ + j * b_dim1;
#line 214 "clagtm.f"
			i__4 = i__ + j * b_dim1;
#line 214 "clagtm.f"
			i__5 = i__ - 1;
#line 214 "clagtm.f"
			i__6 = i__ - 1 + j * x_dim1;
#line 214 "clagtm.f"
			z__4.r = dl[i__5].r * x[i__6].r - dl[i__5].i * x[i__6]
				.i, z__4.i = dl[i__5].r * x[i__6].i + dl[i__5]
				.i * x[i__6].r;
#line 214 "clagtm.f"
			z__3.r = b[i__4].r + z__4.r, z__3.i = b[i__4].i + 
				z__4.i;
#line 214 "clagtm.f"
			i__7 = i__;
#line 214 "clagtm.f"
			i__8 = i__ + j * x_dim1;
#line 214 "clagtm.f"
			z__5.r = d__[i__7].r * x[i__8].r - d__[i__7].i * x[
				i__8].i, z__5.i = d__[i__7].r * x[i__8].i + 
				d__[i__7].i * x[i__8].r;
#line 214 "clagtm.f"
			z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
#line 214 "clagtm.f"
			i__9 = i__;
#line 214 "clagtm.f"
			i__10 = i__ + 1 + j * x_dim1;
#line 214 "clagtm.f"
			z__6.r = du[i__9].r * x[i__10].r - du[i__9].i * x[
				i__10].i, z__6.i = du[i__9].r * x[i__10].i + 
				du[i__9].i * x[i__10].r;
#line 214 "clagtm.f"
			z__1.r = z__2.r + z__6.r, z__1.i = z__2.i + z__6.i;
#line 214 "clagtm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 216 "clagtm.f"
/* L50: */
#line 216 "clagtm.f"
		    }
#line 217 "clagtm.f"
		}
#line 218 "clagtm.f"
/* L60: */
#line 218 "clagtm.f"
	    }
#line 219 "clagtm.f"
	} else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B + A**T * X */

#line 223 "clagtm.f"
	    i__1 = *nrhs;
#line 223 "clagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 224 "clagtm.f"
		if (*n == 1) {
#line 225 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 225 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 225 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 225 "clagtm.f"
		    z__2.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i, 
			    z__2.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4]
			    .r;
#line 225 "clagtm.f"
		    z__1.r = b[i__3].r + z__2.r, z__1.i = b[i__3].i + z__2.i;
#line 225 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 226 "clagtm.f"
		} else {
#line 227 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 227 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 227 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 227 "clagtm.f"
		    z__3.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i, 
			    z__3.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4]
			    .r;
#line 227 "clagtm.f"
		    z__2.r = b[i__3].r + z__3.r, z__2.i = b[i__3].i + z__3.i;
#line 227 "clagtm.f"
		    i__5 = j * x_dim1 + 2;
#line 227 "clagtm.f"
		    z__4.r = dl[1].r * x[i__5].r - dl[1].i * x[i__5].i, 
			    z__4.i = dl[1].r * x[i__5].i + dl[1].i * x[i__5]
			    .r;
#line 227 "clagtm.f"
		    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 227 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 229 "clagtm.f"
		    i__2 = *n + j * b_dim1;
#line 229 "clagtm.f"
		    i__3 = *n + j * b_dim1;
#line 229 "clagtm.f"
		    i__4 = *n - 1;
#line 229 "clagtm.f"
		    i__5 = *n - 1 + j * x_dim1;
#line 229 "clagtm.f"
		    z__3.r = du[i__4].r * x[i__5].r - du[i__4].i * x[i__5].i, 
			    z__3.i = du[i__4].r * x[i__5].i + du[i__4].i * x[
			    i__5].r;
#line 229 "clagtm.f"
		    z__2.r = b[i__3].r + z__3.r, z__2.i = b[i__3].i + z__3.i;
#line 229 "clagtm.f"
		    i__6 = *n;
#line 229 "clagtm.f"
		    i__7 = *n + j * x_dim1;
#line 229 "clagtm.f"
		    z__4.r = d__[i__6].r * x[i__7].r - d__[i__6].i * x[i__7]
			    .i, z__4.i = d__[i__6].r * x[i__7].i + d__[i__6]
			    .i * x[i__7].r;
#line 229 "clagtm.f"
		    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 229 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 231 "clagtm.f"
		    i__2 = *n - 1;
#line 231 "clagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 232 "clagtm.f"
			i__3 = i__ + j * b_dim1;
#line 232 "clagtm.f"
			i__4 = i__ + j * b_dim1;
#line 232 "clagtm.f"
			i__5 = i__ - 1;
#line 232 "clagtm.f"
			i__6 = i__ - 1 + j * x_dim1;
#line 232 "clagtm.f"
			z__4.r = du[i__5].r * x[i__6].r - du[i__5].i * x[i__6]
				.i, z__4.i = du[i__5].r * x[i__6].i + du[i__5]
				.i * x[i__6].r;
#line 232 "clagtm.f"
			z__3.r = b[i__4].r + z__4.r, z__3.i = b[i__4].i + 
				z__4.i;
#line 232 "clagtm.f"
			i__7 = i__;
#line 232 "clagtm.f"
			i__8 = i__ + j * x_dim1;
#line 232 "clagtm.f"
			z__5.r = d__[i__7].r * x[i__8].r - d__[i__7].i * x[
				i__8].i, z__5.i = d__[i__7].r * x[i__8].i + 
				d__[i__7].i * x[i__8].r;
#line 232 "clagtm.f"
			z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
#line 232 "clagtm.f"
			i__9 = i__;
#line 232 "clagtm.f"
			i__10 = i__ + 1 + j * x_dim1;
#line 232 "clagtm.f"
			z__6.r = dl[i__9].r * x[i__10].r - dl[i__9].i * x[
				i__10].i, z__6.i = dl[i__9].r * x[i__10].i + 
				dl[i__9].i * x[i__10].r;
#line 232 "clagtm.f"
			z__1.r = z__2.r + z__6.r, z__1.i = z__2.i + z__6.i;
#line 232 "clagtm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 234 "clagtm.f"
/* L70: */
#line 234 "clagtm.f"
		    }
#line 235 "clagtm.f"
		}
#line 236 "clagtm.f"
/* L80: */
#line 236 "clagtm.f"
	    }
#line 237 "clagtm.f"
	} else if (lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B + A**H * X */

#line 241 "clagtm.f"
	    i__1 = *nrhs;
#line 241 "clagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 242 "clagtm.f"
		if (*n == 1) {
#line 243 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 243 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 243 "clagtm.f"
		    d_cnjg(&z__3, &d__[1]);
#line 243 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 243 "clagtm.f"
		    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, z__2.i =
			     z__3.r * x[i__4].i + z__3.i * x[i__4].r;
#line 243 "clagtm.f"
		    z__1.r = b[i__3].r + z__2.r, z__1.i = b[i__3].i + z__2.i;
#line 243 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 244 "clagtm.f"
		} else {
#line 245 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 245 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 245 "clagtm.f"
		    d_cnjg(&z__4, &d__[1]);
#line 245 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 245 "clagtm.f"
		    z__3.r = z__4.r * x[i__4].r - z__4.i * x[i__4].i, z__3.i =
			     z__4.r * x[i__4].i + z__4.i * x[i__4].r;
#line 245 "clagtm.f"
		    z__2.r = b[i__3].r + z__3.r, z__2.i = b[i__3].i + z__3.i;
#line 245 "clagtm.f"
		    d_cnjg(&z__6, &dl[1]);
#line 245 "clagtm.f"
		    i__5 = j * x_dim1 + 2;
#line 245 "clagtm.f"
		    z__5.r = z__6.r * x[i__5].r - z__6.i * x[i__5].i, z__5.i =
			     z__6.r * x[i__5].i + z__6.i * x[i__5].r;
#line 245 "clagtm.f"
		    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 245 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 247 "clagtm.f"
		    i__2 = *n + j * b_dim1;
#line 247 "clagtm.f"
		    i__3 = *n + j * b_dim1;
#line 247 "clagtm.f"
		    d_cnjg(&z__4, &du[*n - 1]);
#line 247 "clagtm.f"
		    i__4 = *n - 1 + j * x_dim1;
#line 247 "clagtm.f"
		    z__3.r = z__4.r * x[i__4].r - z__4.i * x[i__4].i, z__3.i =
			     z__4.r * x[i__4].i + z__4.i * x[i__4].r;
#line 247 "clagtm.f"
		    z__2.r = b[i__3].r + z__3.r, z__2.i = b[i__3].i + z__3.i;
#line 247 "clagtm.f"
		    d_cnjg(&z__6, &d__[*n]);
#line 247 "clagtm.f"
		    i__5 = *n + j * x_dim1;
#line 247 "clagtm.f"
		    z__5.r = z__6.r * x[i__5].r - z__6.i * x[i__5].i, z__5.i =
			     z__6.r * x[i__5].i + z__6.i * x[i__5].r;
#line 247 "clagtm.f"
		    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 247 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 249 "clagtm.f"
		    i__2 = *n - 1;
#line 249 "clagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 250 "clagtm.f"
			i__3 = i__ + j * b_dim1;
#line 250 "clagtm.f"
			i__4 = i__ + j * b_dim1;
#line 250 "clagtm.f"
			d_cnjg(&z__5, &du[i__ - 1]);
#line 250 "clagtm.f"
			i__5 = i__ - 1 + j * x_dim1;
#line 250 "clagtm.f"
			z__4.r = z__5.r * x[i__5].r - z__5.i * x[i__5].i, 
				z__4.i = z__5.r * x[i__5].i + z__5.i * x[i__5]
				.r;
#line 250 "clagtm.f"
			z__3.r = b[i__4].r + z__4.r, z__3.i = b[i__4].i + 
				z__4.i;
#line 250 "clagtm.f"
			d_cnjg(&z__7, &d__[i__]);
#line 250 "clagtm.f"
			i__6 = i__ + j * x_dim1;
#line 250 "clagtm.f"
			z__6.r = z__7.r * x[i__6].r - z__7.i * x[i__6].i, 
				z__6.i = z__7.r * x[i__6].i + z__7.i * x[i__6]
				.r;
#line 250 "clagtm.f"
			z__2.r = z__3.r + z__6.r, z__2.i = z__3.i + z__6.i;
#line 250 "clagtm.f"
			d_cnjg(&z__9, &dl[i__]);
#line 250 "clagtm.f"
			i__7 = i__ + 1 + j * x_dim1;
#line 250 "clagtm.f"
			z__8.r = z__9.r * x[i__7].r - z__9.i * x[i__7].i, 
				z__8.i = z__9.r * x[i__7].i + z__9.i * x[i__7]
				.r;
#line 250 "clagtm.f"
			z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 250 "clagtm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 254 "clagtm.f"
/* L90: */
#line 254 "clagtm.f"
		    }
#line 255 "clagtm.f"
		}
#line 256 "clagtm.f"
/* L100: */
#line 256 "clagtm.f"
	    }
#line 257 "clagtm.f"
	}
#line 258 "clagtm.f"
    } else if (*alpha == -1.) {
#line 259 "clagtm.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B - A*X */

#line 263 "clagtm.f"
	    i__1 = *nrhs;
#line 263 "clagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 264 "clagtm.f"
		if (*n == 1) {
#line 265 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 265 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 265 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 265 "clagtm.f"
		    z__2.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i, 
			    z__2.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4]
			    .r;
#line 265 "clagtm.f"
		    z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
#line 265 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 266 "clagtm.f"
		} else {
#line 267 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 267 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 267 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 267 "clagtm.f"
		    z__3.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i, 
			    z__3.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4]
			    .r;
#line 267 "clagtm.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 267 "clagtm.f"
		    i__5 = j * x_dim1 + 2;
#line 267 "clagtm.f"
		    z__4.r = du[1].r * x[i__5].r - du[1].i * x[i__5].i, 
			    z__4.i = du[1].r * x[i__5].i + du[1].i * x[i__5]
			    .r;
#line 267 "clagtm.f"
		    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
#line 267 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 269 "clagtm.f"
		    i__2 = *n + j * b_dim1;
#line 269 "clagtm.f"
		    i__3 = *n + j * b_dim1;
#line 269 "clagtm.f"
		    i__4 = *n - 1;
#line 269 "clagtm.f"
		    i__5 = *n - 1 + j * x_dim1;
#line 269 "clagtm.f"
		    z__3.r = dl[i__4].r * x[i__5].r - dl[i__4].i * x[i__5].i, 
			    z__3.i = dl[i__4].r * x[i__5].i + dl[i__4].i * x[
			    i__5].r;
#line 269 "clagtm.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 269 "clagtm.f"
		    i__6 = *n;
#line 269 "clagtm.f"
		    i__7 = *n + j * x_dim1;
#line 269 "clagtm.f"
		    z__4.r = d__[i__6].r * x[i__7].r - d__[i__6].i * x[i__7]
			    .i, z__4.i = d__[i__6].r * x[i__7].i + d__[i__6]
			    .i * x[i__7].r;
#line 269 "clagtm.f"
		    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
#line 269 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 271 "clagtm.f"
		    i__2 = *n - 1;
#line 271 "clagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 272 "clagtm.f"
			i__3 = i__ + j * b_dim1;
#line 272 "clagtm.f"
			i__4 = i__ + j * b_dim1;
#line 272 "clagtm.f"
			i__5 = i__ - 1;
#line 272 "clagtm.f"
			i__6 = i__ - 1 + j * x_dim1;
#line 272 "clagtm.f"
			z__4.r = dl[i__5].r * x[i__6].r - dl[i__5].i * x[i__6]
				.i, z__4.i = dl[i__5].r * x[i__6].i + dl[i__5]
				.i * x[i__6].r;
#line 272 "clagtm.f"
			z__3.r = b[i__4].r - z__4.r, z__3.i = b[i__4].i - 
				z__4.i;
#line 272 "clagtm.f"
			i__7 = i__;
#line 272 "clagtm.f"
			i__8 = i__ + j * x_dim1;
#line 272 "clagtm.f"
			z__5.r = d__[i__7].r * x[i__8].r - d__[i__7].i * x[
				i__8].i, z__5.i = d__[i__7].r * x[i__8].i + 
				d__[i__7].i * x[i__8].r;
#line 272 "clagtm.f"
			z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
#line 272 "clagtm.f"
			i__9 = i__;
#line 272 "clagtm.f"
			i__10 = i__ + 1 + j * x_dim1;
#line 272 "clagtm.f"
			z__6.r = du[i__9].r * x[i__10].r - du[i__9].i * x[
				i__10].i, z__6.i = du[i__9].r * x[i__10].i + 
				du[i__9].i * x[i__10].r;
#line 272 "clagtm.f"
			z__1.r = z__2.r - z__6.r, z__1.i = z__2.i - z__6.i;
#line 272 "clagtm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 274 "clagtm.f"
/* L110: */
#line 274 "clagtm.f"
		    }
#line 275 "clagtm.f"
		}
#line 276 "clagtm.f"
/* L120: */
#line 276 "clagtm.f"
	    }
#line 277 "clagtm.f"
	} else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B - A**T*X */

#line 281 "clagtm.f"
	    i__1 = *nrhs;
#line 281 "clagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 282 "clagtm.f"
		if (*n == 1) {
#line 283 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 283 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 283 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 283 "clagtm.f"
		    z__2.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i, 
			    z__2.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4]
			    .r;
#line 283 "clagtm.f"
		    z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
#line 283 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 284 "clagtm.f"
		} else {
#line 285 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 285 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 285 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 285 "clagtm.f"
		    z__3.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i, 
			    z__3.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4]
			    .r;
#line 285 "clagtm.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 285 "clagtm.f"
		    i__5 = j * x_dim1 + 2;
#line 285 "clagtm.f"
		    z__4.r = dl[1].r * x[i__5].r - dl[1].i * x[i__5].i, 
			    z__4.i = dl[1].r * x[i__5].i + dl[1].i * x[i__5]
			    .r;
#line 285 "clagtm.f"
		    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
#line 285 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 287 "clagtm.f"
		    i__2 = *n + j * b_dim1;
#line 287 "clagtm.f"
		    i__3 = *n + j * b_dim1;
#line 287 "clagtm.f"
		    i__4 = *n - 1;
#line 287 "clagtm.f"
		    i__5 = *n - 1 + j * x_dim1;
#line 287 "clagtm.f"
		    z__3.r = du[i__4].r * x[i__5].r - du[i__4].i * x[i__5].i, 
			    z__3.i = du[i__4].r * x[i__5].i + du[i__4].i * x[
			    i__5].r;
#line 287 "clagtm.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 287 "clagtm.f"
		    i__6 = *n;
#line 287 "clagtm.f"
		    i__7 = *n + j * x_dim1;
#line 287 "clagtm.f"
		    z__4.r = d__[i__6].r * x[i__7].r - d__[i__6].i * x[i__7]
			    .i, z__4.i = d__[i__6].r * x[i__7].i + d__[i__6]
			    .i * x[i__7].r;
#line 287 "clagtm.f"
		    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
#line 287 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 289 "clagtm.f"
		    i__2 = *n - 1;
#line 289 "clagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 290 "clagtm.f"
			i__3 = i__ + j * b_dim1;
#line 290 "clagtm.f"
			i__4 = i__ + j * b_dim1;
#line 290 "clagtm.f"
			i__5 = i__ - 1;
#line 290 "clagtm.f"
			i__6 = i__ - 1 + j * x_dim1;
#line 290 "clagtm.f"
			z__4.r = du[i__5].r * x[i__6].r - du[i__5].i * x[i__6]
				.i, z__4.i = du[i__5].r * x[i__6].i + du[i__5]
				.i * x[i__6].r;
#line 290 "clagtm.f"
			z__3.r = b[i__4].r - z__4.r, z__3.i = b[i__4].i - 
				z__4.i;
#line 290 "clagtm.f"
			i__7 = i__;
#line 290 "clagtm.f"
			i__8 = i__ + j * x_dim1;
#line 290 "clagtm.f"
			z__5.r = d__[i__7].r * x[i__8].r - d__[i__7].i * x[
				i__8].i, z__5.i = d__[i__7].r * x[i__8].i + 
				d__[i__7].i * x[i__8].r;
#line 290 "clagtm.f"
			z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
#line 290 "clagtm.f"
			i__9 = i__;
#line 290 "clagtm.f"
			i__10 = i__ + 1 + j * x_dim1;
#line 290 "clagtm.f"
			z__6.r = dl[i__9].r * x[i__10].r - dl[i__9].i * x[
				i__10].i, z__6.i = dl[i__9].r * x[i__10].i + 
				dl[i__9].i * x[i__10].r;
#line 290 "clagtm.f"
			z__1.r = z__2.r - z__6.r, z__1.i = z__2.i - z__6.i;
#line 290 "clagtm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 292 "clagtm.f"
/* L130: */
#line 292 "clagtm.f"
		    }
#line 293 "clagtm.f"
		}
#line 294 "clagtm.f"
/* L140: */
#line 294 "clagtm.f"
	    }
#line 295 "clagtm.f"
	} else if (lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B - A**H*X */

#line 299 "clagtm.f"
	    i__1 = *nrhs;
#line 299 "clagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 300 "clagtm.f"
		if (*n == 1) {
#line 301 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 301 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 301 "clagtm.f"
		    d_cnjg(&z__3, &d__[1]);
#line 301 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 301 "clagtm.f"
		    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, z__2.i =
			     z__3.r * x[i__4].i + z__3.i * x[i__4].r;
#line 301 "clagtm.f"
		    z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
#line 301 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 302 "clagtm.f"
		} else {
#line 303 "clagtm.f"
		    i__2 = j * b_dim1 + 1;
#line 303 "clagtm.f"
		    i__3 = j * b_dim1 + 1;
#line 303 "clagtm.f"
		    d_cnjg(&z__4, &d__[1]);
#line 303 "clagtm.f"
		    i__4 = j * x_dim1 + 1;
#line 303 "clagtm.f"
		    z__3.r = z__4.r * x[i__4].r - z__4.i * x[i__4].i, z__3.i =
			     z__4.r * x[i__4].i + z__4.i * x[i__4].r;
#line 303 "clagtm.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 303 "clagtm.f"
		    d_cnjg(&z__6, &dl[1]);
#line 303 "clagtm.f"
		    i__5 = j * x_dim1 + 2;
#line 303 "clagtm.f"
		    z__5.r = z__6.r * x[i__5].r - z__6.i * x[i__5].i, z__5.i =
			     z__6.r * x[i__5].i + z__6.i * x[i__5].r;
#line 303 "clagtm.f"
		    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - z__5.i;
#line 303 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 305 "clagtm.f"
		    i__2 = *n + j * b_dim1;
#line 305 "clagtm.f"
		    i__3 = *n + j * b_dim1;
#line 305 "clagtm.f"
		    d_cnjg(&z__4, &du[*n - 1]);
#line 305 "clagtm.f"
		    i__4 = *n - 1 + j * x_dim1;
#line 305 "clagtm.f"
		    z__3.r = z__4.r * x[i__4].r - z__4.i * x[i__4].i, z__3.i =
			     z__4.r * x[i__4].i + z__4.i * x[i__4].r;
#line 305 "clagtm.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 305 "clagtm.f"
		    d_cnjg(&z__6, &d__[*n]);
#line 305 "clagtm.f"
		    i__5 = *n + j * x_dim1;
#line 305 "clagtm.f"
		    z__5.r = z__6.r * x[i__5].r - z__6.i * x[i__5].i, z__5.i =
			     z__6.r * x[i__5].i + z__6.i * x[i__5].r;
#line 305 "clagtm.f"
		    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - z__5.i;
#line 305 "clagtm.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 307 "clagtm.f"
		    i__2 = *n - 1;
#line 307 "clagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 308 "clagtm.f"
			i__3 = i__ + j * b_dim1;
#line 308 "clagtm.f"
			i__4 = i__ + j * b_dim1;
#line 308 "clagtm.f"
			d_cnjg(&z__5, &du[i__ - 1]);
#line 308 "clagtm.f"
			i__5 = i__ - 1 + j * x_dim1;
#line 308 "clagtm.f"
			z__4.r = z__5.r * x[i__5].r - z__5.i * x[i__5].i, 
				z__4.i = z__5.r * x[i__5].i + z__5.i * x[i__5]
				.r;
#line 308 "clagtm.f"
			z__3.r = b[i__4].r - z__4.r, z__3.i = b[i__4].i - 
				z__4.i;
#line 308 "clagtm.f"
			d_cnjg(&z__7, &d__[i__]);
#line 308 "clagtm.f"
			i__6 = i__ + j * x_dim1;
#line 308 "clagtm.f"
			z__6.r = z__7.r * x[i__6].r - z__7.i * x[i__6].i, 
				z__6.i = z__7.r * x[i__6].i + z__7.i * x[i__6]
				.r;
#line 308 "clagtm.f"
			z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 308 "clagtm.f"
			d_cnjg(&z__9, &dl[i__]);
#line 308 "clagtm.f"
			i__7 = i__ + 1 + j * x_dim1;
#line 308 "clagtm.f"
			z__8.r = z__9.r * x[i__7].r - z__9.i * x[i__7].i, 
				z__8.i = z__9.r * x[i__7].i + z__9.i * x[i__7]
				.r;
#line 308 "clagtm.f"
			z__1.r = z__2.r - z__8.r, z__1.i = z__2.i - z__8.i;
#line 308 "clagtm.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 312 "clagtm.f"
/* L150: */
#line 312 "clagtm.f"
		    }
#line 313 "clagtm.f"
		}
#line 314 "clagtm.f"
/* L160: */
#line 314 "clagtm.f"
	    }
#line 315 "clagtm.f"
	}
#line 316 "clagtm.f"
    }
#line 317 "clagtm.f"
    return 0;

/*     End of CLAGTM */

} /* clagtm_ */

