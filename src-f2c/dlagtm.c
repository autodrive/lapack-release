#line 1 "dlagtm.f"
/* dlagtm.f -- translated by f2c (version 20100827).
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

#line 1 "dlagtm.f"
/* > \brief \b DLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matr
ix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAGTM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagtm.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagtm.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagtm.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, */
/*                          B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            LDB, LDX, N, NRHS */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAGTM performs a matrix-vector product of the form */
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
/* >          = 'T':  Transpose,    B := alpha * A'* X + beta * B */
/* >          = 'C':  Conjugate transpose = Transpose */
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
/* >          ALPHA is DOUBLE PRECISION */
/* >          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise, */
/* >          it is assumed to be 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) sub-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) super-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
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
/* >          BETA is DOUBLE PRECISION */
/* >          The scalar beta.  BETA must be 0., 1., or -1.; otherwise, */
/* >          it is assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlagtm_(char *trans, integer *n, integer *nrhs, 
	doublereal *alpha, doublereal *dl, doublereal *d__, doublereal *du, 
	doublereal *x, integer *ldx, doublereal *beta, doublereal *b, integer 
	*ldb, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;

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
/*     .. Executable Statements .. */

#line 178 "dlagtm.f"
    /* Parameter adjustments */
#line 178 "dlagtm.f"
    --dl;
#line 178 "dlagtm.f"
    --d__;
#line 178 "dlagtm.f"
    --du;
#line 178 "dlagtm.f"
    x_dim1 = *ldx;
#line 178 "dlagtm.f"
    x_offset = 1 + x_dim1;
#line 178 "dlagtm.f"
    x -= x_offset;
#line 178 "dlagtm.f"
    b_dim1 = *ldb;
#line 178 "dlagtm.f"
    b_offset = 1 + b_dim1;
#line 178 "dlagtm.f"
    b -= b_offset;
#line 178 "dlagtm.f"

#line 178 "dlagtm.f"
    /* Function Body */
#line 178 "dlagtm.f"
    if (*n == 0) {
#line 178 "dlagtm.f"
	return 0;
#line 178 "dlagtm.f"
    }

/*     Multiply B by BETA if BETA.NE.1. */

#line 183 "dlagtm.f"
    if (*beta == 0.) {
#line 184 "dlagtm.f"
	i__1 = *nrhs;
#line 184 "dlagtm.f"
	for (j = 1; j <= i__1; ++j) {
#line 185 "dlagtm.f"
	    i__2 = *n;
#line 185 "dlagtm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 186 "dlagtm.f"
		b[i__ + j * b_dim1] = 0.;
#line 187 "dlagtm.f"
/* L10: */
#line 187 "dlagtm.f"
	    }
#line 188 "dlagtm.f"
/* L20: */
#line 188 "dlagtm.f"
	}
#line 189 "dlagtm.f"
    } else if (*beta == -1.) {
#line 190 "dlagtm.f"
	i__1 = *nrhs;
#line 190 "dlagtm.f"
	for (j = 1; j <= i__1; ++j) {
#line 191 "dlagtm.f"
	    i__2 = *n;
#line 191 "dlagtm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 192 "dlagtm.f"
		b[i__ + j * b_dim1] = -b[i__ + j * b_dim1];
#line 193 "dlagtm.f"
/* L30: */
#line 193 "dlagtm.f"
	    }
#line 194 "dlagtm.f"
/* L40: */
#line 194 "dlagtm.f"
	}
#line 195 "dlagtm.f"
    }

#line 197 "dlagtm.f"
    if (*alpha == 1.) {
#line 198 "dlagtm.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B + A*X */

#line 202 "dlagtm.f"
	    i__1 = *nrhs;
#line 202 "dlagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 203 "dlagtm.f"
		if (*n == 1) {
#line 204 "dlagtm.f"
		    b[j * b_dim1 + 1] += d__[1] * x[j * x_dim1 + 1];
#line 205 "dlagtm.f"
		} else {
#line 206 "dlagtm.f"
		    b[j * b_dim1 + 1] = b[j * b_dim1 + 1] + d__[1] * x[j * 
			    x_dim1 + 1] + du[1] * x[j * x_dim1 + 2];
#line 208 "dlagtm.f"
		    b[*n + j * b_dim1] = b[*n + j * b_dim1] + dl[*n - 1] * x[*
			    n - 1 + j * x_dim1] + d__[*n] * x[*n + j * x_dim1]
			    ;
#line 210 "dlagtm.f"
		    i__2 = *n - 1;
#line 210 "dlagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 211 "dlagtm.f"
			b[i__ + j * b_dim1] = b[i__ + j * b_dim1] + dl[i__ - 
				1] * x[i__ - 1 + j * x_dim1] + d__[i__] * x[
				i__ + j * x_dim1] + du[i__] * x[i__ + 1 + j * 
				x_dim1];
#line 213 "dlagtm.f"
/* L50: */
#line 213 "dlagtm.f"
		    }
#line 214 "dlagtm.f"
		}
#line 215 "dlagtm.f"
/* L60: */
#line 215 "dlagtm.f"
	    }
#line 216 "dlagtm.f"
	} else {

/*           Compute B := B + A**T*X */

#line 220 "dlagtm.f"
	    i__1 = *nrhs;
#line 220 "dlagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 221 "dlagtm.f"
		if (*n == 1) {
#line 222 "dlagtm.f"
		    b[j * b_dim1 + 1] += d__[1] * x[j * x_dim1 + 1];
#line 223 "dlagtm.f"
		} else {
#line 224 "dlagtm.f"
		    b[j * b_dim1 + 1] = b[j * b_dim1 + 1] + d__[1] * x[j * 
			    x_dim1 + 1] + dl[1] * x[j * x_dim1 + 2];
#line 226 "dlagtm.f"
		    b[*n + j * b_dim1] = b[*n + j * b_dim1] + du[*n - 1] * x[*
			    n - 1 + j * x_dim1] + d__[*n] * x[*n + j * x_dim1]
			    ;
#line 228 "dlagtm.f"
		    i__2 = *n - 1;
#line 228 "dlagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 229 "dlagtm.f"
			b[i__ + j * b_dim1] = b[i__ + j * b_dim1] + du[i__ - 
				1] * x[i__ - 1 + j * x_dim1] + d__[i__] * x[
				i__ + j * x_dim1] + dl[i__] * x[i__ + 1 + j * 
				x_dim1];
#line 231 "dlagtm.f"
/* L70: */
#line 231 "dlagtm.f"
		    }
#line 232 "dlagtm.f"
		}
#line 233 "dlagtm.f"
/* L80: */
#line 233 "dlagtm.f"
	    }
#line 234 "dlagtm.f"
	}
#line 235 "dlagtm.f"
    } else if (*alpha == -1.) {
#line 236 "dlagtm.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B - A*X */

#line 240 "dlagtm.f"
	    i__1 = *nrhs;
#line 240 "dlagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 241 "dlagtm.f"
		if (*n == 1) {
#line 242 "dlagtm.f"
		    b[j * b_dim1 + 1] -= d__[1] * x[j * x_dim1 + 1];
#line 243 "dlagtm.f"
		} else {
#line 244 "dlagtm.f"
		    b[j * b_dim1 + 1] = b[j * b_dim1 + 1] - d__[1] * x[j * 
			    x_dim1 + 1] - du[1] * x[j * x_dim1 + 2];
#line 246 "dlagtm.f"
		    b[*n + j * b_dim1] = b[*n + j * b_dim1] - dl[*n - 1] * x[*
			    n - 1 + j * x_dim1] - d__[*n] * x[*n + j * x_dim1]
			    ;
#line 248 "dlagtm.f"
		    i__2 = *n - 1;
#line 248 "dlagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 249 "dlagtm.f"
			b[i__ + j * b_dim1] = b[i__ + j * b_dim1] - dl[i__ - 
				1] * x[i__ - 1 + j * x_dim1] - d__[i__] * x[
				i__ + j * x_dim1] - du[i__] * x[i__ + 1 + j * 
				x_dim1];
#line 251 "dlagtm.f"
/* L90: */
#line 251 "dlagtm.f"
		    }
#line 252 "dlagtm.f"
		}
#line 253 "dlagtm.f"
/* L100: */
#line 253 "dlagtm.f"
	    }
#line 254 "dlagtm.f"
	} else {

/*           Compute B := B - A**T*X */

#line 258 "dlagtm.f"
	    i__1 = *nrhs;
#line 258 "dlagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 259 "dlagtm.f"
		if (*n == 1) {
#line 260 "dlagtm.f"
		    b[j * b_dim1 + 1] -= d__[1] * x[j * x_dim1 + 1];
#line 261 "dlagtm.f"
		} else {
#line 262 "dlagtm.f"
		    b[j * b_dim1 + 1] = b[j * b_dim1 + 1] - d__[1] * x[j * 
			    x_dim1 + 1] - dl[1] * x[j * x_dim1 + 2];
#line 264 "dlagtm.f"
		    b[*n + j * b_dim1] = b[*n + j * b_dim1] - du[*n - 1] * x[*
			    n - 1 + j * x_dim1] - d__[*n] * x[*n + j * x_dim1]
			    ;
#line 266 "dlagtm.f"
		    i__2 = *n - 1;
#line 266 "dlagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 267 "dlagtm.f"
			b[i__ + j * b_dim1] = b[i__ + j * b_dim1] - du[i__ - 
				1] * x[i__ - 1 + j * x_dim1] - d__[i__] * x[
				i__ + j * x_dim1] - dl[i__] * x[i__ + 1 + j * 
				x_dim1];
#line 269 "dlagtm.f"
/* L110: */
#line 269 "dlagtm.f"
		    }
#line 270 "dlagtm.f"
		}
#line 271 "dlagtm.f"
/* L120: */
#line 271 "dlagtm.f"
	    }
#line 272 "dlagtm.f"
	}
#line 273 "dlagtm.f"
    }
#line 274 "dlagtm.f"
    return 0;

/*     End of DLAGTM */

} /* dlagtm_ */

