#line 1 "slagtm.f"
/* slagtm.f -- translated by f2c (version 20100827).
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

#line 1 "slagtm.f"
/* > \brief \b SLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matr
ix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAGTM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagtm.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagtm.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagtm.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, */
/*                          B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            LDB, LDX, N, NRHS */
/*       REAL               ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAGTM performs a matrix-vector product of the form */
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
/* >          ALPHA is REAL */
/* >          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise, */
/* >          it is assumed to be 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is REAL array, dimension (N-1) */
/* >          The (n-1) sub-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is REAL array, dimension (N-1) */
/* >          The (n-1) super-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array, dimension (LDX,NRHS) */
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
/* >          B is REAL array, dimension (LDB,NRHS) */
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

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slagtm_(char *trans, integer *n, integer *nrhs, 
	doublereal *alpha, doublereal *dl, doublereal *d__, doublereal *du, 
	doublereal *x, integer *ldx, doublereal *beta, doublereal *b, integer 
	*ldb, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 178 "slagtm.f"
    /* Parameter adjustments */
#line 178 "slagtm.f"
    --dl;
#line 178 "slagtm.f"
    --d__;
#line 178 "slagtm.f"
    --du;
#line 178 "slagtm.f"
    x_dim1 = *ldx;
#line 178 "slagtm.f"
    x_offset = 1 + x_dim1;
#line 178 "slagtm.f"
    x -= x_offset;
#line 178 "slagtm.f"
    b_dim1 = *ldb;
#line 178 "slagtm.f"
    b_offset = 1 + b_dim1;
#line 178 "slagtm.f"
    b -= b_offset;
#line 178 "slagtm.f"

#line 178 "slagtm.f"
    /* Function Body */
#line 178 "slagtm.f"
    if (*n == 0) {
#line 178 "slagtm.f"
	return 0;
#line 178 "slagtm.f"
    }

/*     Multiply B by BETA if BETA.NE.1. */

#line 183 "slagtm.f"
    if (*beta == 0.) {
#line 184 "slagtm.f"
	i__1 = *nrhs;
#line 184 "slagtm.f"
	for (j = 1; j <= i__1; ++j) {
#line 185 "slagtm.f"
	    i__2 = *n;
#line 185 "slagtm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 186 "slagtm.f"
		b[i__ + j * b_dim1] = 0.;
#line 187 "slagtm.f"
/* L10: */
#line 187 "slagtm.f"
	    }
#line 188 "slagtm.f"
/* L20: */
#line 188 "slagtm.f"
	}
#line 189 "slagtm.f"
    } else if (*beta == -1.) {
#line 190 "slagtm.f"
	i__1 = *nrhs;
#line 190 "slagtm.f"
	for (j = 1; j <= i__1; ++j) {
#line 191 "slagtm.f"
	    i__2 = *n;
#line 191 "slagtm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 192 "slagtm.f"
		b[i__ + j * b_dim1] = -b[i__ + j * b_dim1];
#line 193 "slagtm.f"
/* L30: */
#line 193 "slagtm.f"
	    }
#line 194 "slagtm.f"
/* L40: */
#line 194 "slagtm.f"
	}
#line 195 "slagtm.f"
    }

#line 197 "slagtm.f"
    if (*alpha == 1.) {
#line 198 "slagtm.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B + A*X */

#line 202 "slagtm.f"
	    i__1 = *nrhs;
#line 202 "slagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 203 "slagtm.f"
		if (*n == 1) {
#line 204 "slagtm.f"
		    b[j * b_dim1 + 1] += d__[1] * x[j * x_dim1 + 1];
#line 205 "slagtm.f"
		} else {
#line 206 "slagtm.f"
		    b[j * b_dim1 + 1] = b[j * b_dim1 + 1] + d__[1] * x[j * 
			    x_dim1 + 1] + du[1] * x[j * x_dim1 + 2];
#line 208 "slagtm.f"
		    b[*n + j * b_dim1] = b[*n + j * b_dim1] + dl[*n - 1] * x[*
			    n - 1 + j * x_dim1] + d__[*n] * x[*n + j * x_dim1]
			    ;
#line 210 "slagtm.f"
		    i__2 = *n - 1;
#line 210 "slagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 211 "slagtm.f"
			b[i__ + j * b_dim1] = b[i__ + j * b_dim1] + dl[i__ - 
				1] * x[i__ - 1 + j * x_dim1] + d__[i__] * x[
				i__ + j * x_dim1] + du[i__] * x[i__ + 1 + j * 
				x_dim1];
#line 213 "slagtm.f"
/* L50: */
#line 213 "slagtm.f"
		    }
#line 214 "slagtm.f"
		}
#line 215 "slagtm.f"
/* L60: */
#line 215 "slagtm.f"
	    }
#line 216 "slagtm.f"
	} else {

/*           Compute B := B + A**T*X */

#line 220 "slagtm.f"
	    i__1 = *nrhs;
#line 220 "slagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 221 "slagtm.f"
		if (*n == 1) {
#line 222 "slagtm.f"
		    b[j * b_dim1 + 1] += d__[1] * x[j * x_dim1 + 1];
#line 223 "slagtm.f"
		} else {
#line 224 "slagtm.f"
		    b[j * b_dim1 + 1] = b[j * b_dim1 + 1] + d__[1] * x[j * 
			    x_dim1 + 1] + dl[1] * x[j * x_dim1 + 2];
#line 226 "slagtm.f"
		    b[*n + j * b_dim1] = b[*n + j * b_dim1] + du[*n - 1] * x[*
			    n - 1 + j * x_dim1] + d__[*n] * x[*n + j * x_dim1]
			    ;
#line 228 "slagtm.f"
		    i__2 = *n - 1;
#line 228 "slagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 229 "slagtm.f"
			b[i__ + j * b_dim1] = b[i__ + j * b_dim1] + du[i__ - 
				1] * x[i__ - 1 + j * x_dim1] + d__[i__] * x[
				i__ + j * x_dim1] + dl[i__] * x[i__ + 1 + j * 
				x_dim1];
#line 231 "slagtm.f"
/* L70: */
#line 231 "slagtm.f"
		    }
#line 232 "slagtm.f"
		}
#line 233 "slagtm.f"
/* L80: */
#line 233 "slagtm.f"
	    }
#line 234 "slagtm.f"
	}
#line 235 "slagtm.f"
    } else if (*alpha == -1.) {
#line 236 "slagtm.f"
	if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*           Compute B := B - A*X */

#line 240 "slagtm.f"
	    i__1 = *nrhs;
#line 240 "slagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 241 "slagtm.f"
		if (*n == 1) {
#line 242 "slagtm.f"
		    b[j * b_dim1 + 1] -= d__[1] * x[j * x_dim1 + 1];
#line 243 "slagtm.f"
		} else {
#line 244 "slagtm.f"
		    b[j * b_dim1 + 1] = b[j * b_dim1 + 1] - d__[1] * x[j * 
			    x_dim1 + 1] - du[1] * x[j * x_dim1 + 2];
#line 246 "slagtm.f"
		    b[*n + j * b_dim1] = b[*n + j * b_dim1] - dl[*n - 1] * x[*
			    n - 1 + j * x_dim1] - d__[*n] * x[*n + j * x_dim1]
			    ;
#line 248 "slagtm.f"
		    i__2 = *n - 1;
#line 248 "slagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 249 "slagtm.f"
			b[i__ + j * b_dim1] = b[i__ + j * b_dim1] - dl[i__ - 
				1] * x[i__ - 1 + j * x_dim1] - d__[i__] * x[
				i__ + j * x_dim1] - du[i__] * x[i__ + 1 + j * 
				x_dim1];
#line 251 "slagtm.f"
/* L90: */
#line 251 "slagtm.f"
		    }
#line 252 "slagtm.f"
		}
#line 253 "slagtm.f"
/* L100: */
#line 253 "slagtm.f"
	    }
#line 254 "slagtm.f"
	} else {

/*           Compute B := B - A**T*X */

#line 258 "slagtm.f"
	    i__1 = *nrhs;
#line 258 "slagtm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 259 "slagtm.f"
		if (*n == 1) {
#line 260 "slagtm.f"
		    b[j * b_dim1 + 1] -= d__[1] * x[j * x_dim1 + 1];
#line 261 "slagtm.f"
		} else {
#line 262 "slagtm.f"
		    b[j * b_dim1 + 1] = b[j * b_dim1 + 1] - d__[1] * x[j * 
			    x_dim1 + 1] - dl[1] * x[j * x_dim1 + 2];
#line 264 "slagtm.f"
		    b[*n + j * b_dim1] = b[*n + j * b_dim1] - du[*n - 1] * x[*
			    n - 1 + j * x_dim1] - d__[*n] * x[*n + j * x_dim1]
			    ;
#line 266 "slagtm.f"
		    i__2 = *n - 1;
#line 266 "slagtm.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 267 "slagtm.f"
			b[i__ + j * b_dim1] = b[i__ + j * b_dim1] - du[i__ - 
				1] * x[i__ - 1 + j * x_dim1] - d__[i__] * x[
				i__ + j * x_dim1] - dl[i__] * x[i__ + 1 + j * 
				x_dim1];
#line 269 "slagtm.f"
/* L110: */
#line 269 "slagtm.f"
		    }
#line 270 "slagtm.f"
		}
#line 271 "slagtm.f"
/* L120: */
#line 271 "slagtm.f"
	    }
#line 272 "slagtm.f"
	}
#line 273 "slagtm.f"
    }
#line 274 "slagtm.f"
    return 0;

/*     End of SLAGTM */

} /* slagtm_ */

