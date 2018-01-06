#line 1 "sgtts2.f"
/* sgtts2.f -- translated by f2c (version 20100827).
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

#line 1 "sgtts2.f"
/* > \brief \b SGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization
 computed by sgttrf. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGTTS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtts2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtts2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtts2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ITRANS, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGTTS2 solves one of the systems of equations */
/* >    A*X = B  or  A**T*X = B, */
/* > with a tridiagonal matrix A using the LU factorization computed */
/* > by SGTTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITRANS */
/* > \verbatim */
/* >          ITRANS is INTEGER */
/* >          Specifies the form of the system of equations. */
/* >          = 0:  A * X = B  (No transpose) */
/* >          = 1:  A**T* X = B  (Transpose) */
/* >          = 2:  A**T* X = B  (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is REAL array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is REAL array, dimension (N-1) */
/* >          The (n-1) elements of the first super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is REAL array, dimension (N-2) */
/* >          The (n-2) elements of the second super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= n, row i of the matrix was */
/* >          interchanged with row IPIV(i).  IPIV(i) will always be either */
/* >          i or i+1; IPIV(i) = i indicates a row interchange was not */
/* >          required. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
/* >          On entry, the matrix of right hand side vectors B. */
/* >          On exit, B is overwritten by the solution vectors X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgtts2_(integer *itrans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, 
	integer *ipiv, doublereal *b, integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ip;
    static doublereal temp;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 154 "sgtts2.f"
    /* Parameter adjustments */
#line 154 "sgtts2.f"
    --dl;
#line 154 "sgtts2.f"
    --d__;
#line 154 "sgtts2.f"
    --du;
#line 154 "sgtts2.f"
    --du2;
#line 154 "sgtts2.f"
    --ipiv;
#line 154 "sgtts2.f"
    b_dim1 = *ldb;
#line 154 "sgtts2.f"
    b_offset = 1 + b_dim1;
#line 154 "sgtts2.f"
    b -= b_offset;
#line 154 "sgtts2.f"

#line 154 "sgtts2.f"
    /* Function Body */
#line 154 "sgtts2.f"
    if (*n == 0 || *nrhs == 0) {
#line 154 "sgtts2.f"
	return 0;
#line 154 "sgtts2.f"
    }

#line 157 "sgtts2.f"
    if (*itrans == 0) {

/*        Solve A*X = B using the LU factorization of A, */
/*        overwriting each right hand side vector with its solution. */

#line 162 "sgtts2.f"
	if (*nrhs <= 1) {
#line 163 "sgtts2.f"
	    j = 1;
#line 164 "sgtts2.f"
L10:

/*           Solve L*x = b. */

#line 168 "sgtts2.f"
	    i__1 = *n - 1;
#line 168 "sgtts2.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 169 "sgtts2.f"
		ip = ipiv[i__];
#line 170 "sgtts2.f"
		temp = b[i__ + 1 - ip + i__ + j * b_dim1] - dl[i__] * b[ip + 
			j * b_dim1];
#line 171 "sgtts2.f"
		b[i__ + j * b_dim1] = b[ip + j * b_dim1];
#line 172 "sgtts2.f"
		b[i__ + 1 + j * b_dim1] = temp;
#line 173 "sgtts2.f"
/* L20: */
#line 173 "sgtts2.f"
	    }

/*           Solve U*x = b. */

#line 177 "sgtts2.f"
	    b[*n + j * b_dim1] /= d__[*n];
#line 178 "sgtts2.f"
	    if (*n > 1) {
#line 178 "sgtts2.f"
		b[*n - 1 + j * b_dim1] = (b[*n - 1 + j * b_dim1] - du[*n - 1] 
			* b[*n + j * b_dim1]) / d__[*n - 1];
#line 178 "sgtts2.f"
	    }
#line 181 "sgtts2.f"
	    for (i__ = *n - 2; i__ >= 1; --i__) {
#line 182 "sgtts2.f"
		b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ 
			+ 1 + j * b_dim1] - du2[i__] * b[i__ + 2 + j * b_dim1]
			) / d__[i__];
#line 184 "sgtts2.f"
/* L30: */
#line 184 "sgtts2.f"
	    }
#line 185 "sgtts2.f"
	    if (j < *nrhs) {
#line 186 "sgtts2.f"
		++j;
#line 187 "sgtts2.f"
		goto L10;
#line 188 "sgtts2.f"
	    }
#line 189 "sgtts2.f"
	} else {
#line 190 "sgtts2.f"
	    i__1 = *nrhs;
#line 190 "sgtts2.f"
	    for (j = 1; j <= i__1; ++j) {

/*              Solve L*x = b. */

#line 194 "sgtts2.f"
		i__2 = *n - 1;
#line 194 "sgtts2.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 195 "sgtts2.f"
		    if (ipiv[i__] == i__) {
#line 196 "sgtts2.f"
			b[i__ + 1 + j * b_dim1] -= dl[i__] * b[i__ + j * 
				b_dim1];
#line 197 "sgtts2.f"
		    } else {
#line 198 "sgtts2.f"
			temp = b[i__ + j * b_dim1];
#line 199 "sgtts2.f"
			b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
#line 200 "sgtts2.f"
			b[i__ + 1 + j * b_dim1] = temp - dl[i__] * b[i__ + j *
				 b_dim1];
#line 201 "sgtts2.f"
		    }
#line 202 "sgtts2.f"
/* L40: */
#line 202 "sgtts2.f"
		}

/*              Solve U*x = b. */

#line 206 "sgtts2.f"
		b[*n + j * b_dim1] /= d__[*n];
#line 207 "sgtts2.f"
		if (*n > 1) {
#line 207 "sgtts2.f"
		    b[*n - 1 + j * b_dim1] = (b[*n - 1 + j * b_dim1] - du[*n 
			    - 1] * b[*n + j * b_dim1]) / d__[*n - 1];
#line 207 "sgtts2.f"
		}
#line 210 "sgtts2.f"
		for (i__ = *n - 2; i__ >= 1; --i__) {
#line 211 "sgtts2.f"
		    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[
			    i__ + 1 + j * b_dim1] - du2[i__] * b[i__ + 2 + j *
			     b_dim1]) / d__[i__];
#line 213 "sgtts2.f"
/* L50: */
#line 213 "sgtts2.f"
		}
#line 214 "sgtts2.f"
/* L60: */
#line 214 "sgtts2.f"
	    }
#line 215 "sgtts2.f"
	}
#line 216 "sgtts2.f"
    } else {

/*        Solve A**T * X = B. */

#line 220 "sgtts2.f"
	if (*nrhs <= 1) {

/*           Solve U**T*x = b. */

#line 224 "sgtts2.f"
	    j = 1;
#line 225 "sgtts2.f"
L70:
#line 226 "sgtts2.f"
	    b[j * b_dim1 + 1] /= d__[1];
#line 227 "sgtts2.f"
	    if (*n > 1) {
#line 227 "sgtts2.f"
		b[j * b_dim1 + 2] = (b[j * b_dim1 + 2] - du[1] * b[j * b_dim1 
			+ 1]) / d__[2];
#line 227 "sgtts2.f"
	    }
#line 229 "sgtts2.f"
	    i__1 = *n;
#line 229 "sgtts2.f"
	    for (i__ = 3; i__ <= i__1; ++i__) {
#line 230 "sgtts2.f"
		b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__ - 1] * b[
			i__ - 1 + j * b_dim1] - du2[i__ - 2] * b[i__ - 2 + j *
			 b_dim1]) / d__[i__];
#line 232 "sgtts2.f"
/* L80: */
#line 232 "sgtts2.f"
	    }

/*           Solve L**T*x = b. */

#line 236 "sgtts2.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 237 "sgtts2.f"
		ip = ipiv[i__];
#line 238 "sgtts2.f"
		temp = b[i__ + j * b_dim1] - dl[i__] * b[i__ + 1 + j * b_dim1]
			;
#line 239 "sgtts2.f"
		b[i__ + j * b_dim1] = b[ip + j * b_dim1];
#line 240 "sgtts2.f"
		b[ip + j * b_dim1] = temp;
#line 241 "sgtts2.f"
/* L90: */
#line 241 "sgtts2.f"
	    }
#line 242 "sgtts2.f"
	    if (j < *nrhs) {
#line 243 "sgtts2.f"
		++j;
#line 244 "sgtts2.f"
		goto L70;
#line 245 "sgtts2.f"
	    }

#line 247 "sgtts2.f"
	} else {
#line 248 "sgtts2.f"
	    i__1 = *nrhs;
#line 248 "sgtts2.f"
	    for (j = 1; j <= i__1; ++j) {

/*              Solve U**T*x = b. */

#line 252 "sgtts2.f"
		b[j * b_dim1 + 1] /= d__[1];
#line 253 "sgtts2.f"
		if (*n > 1) {
#line 253 "sgtts2.f"
		    b[j * b_dim1 + 2] = (b[j * b_dim1 + 2] - du[1] * b[j * 
			    b_dim1 + 1]) / d__[2];
#line 253 "sgtts2.f"
		}
#line 255 "sgtts2.f"
		i__2 = *n;
#line 255 "sgtts2.f"
		for (i__ = 3; i__ <= i__2; ++i__) {
#line 256 "sgtts2.f"
		    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__ - 1] *
			     b[i__ - 1 + j * b_dim1] - du2[i__ - 2] * b[i__ - 
			    2 + j * b_dim1]) / d__[i__];
#line 258 "sgtts2.f"
/* L100: */
#line 258 "sgtts2.f"
		}
#line 259 "sgtts2.f"
		for (i__ = *n - 1; i__ >= 1; --i__) {
#line 260 "sgtts2.f"
		    if (ipiv[i__] == i__) {
#line 261 "sgtts2.f"
			b[i__ + j * b_dim1] -= dl[i__] * b[i__ + 1 + j * 
				b_dim1];
#line 262 "sgtts2.f"
		    } else {
#line 263 "sgtts2.f"
			temp = b[i__ + 1 + j * b_dim1];
#line 264 "sgtts2.f"
			b[i__ + 1 + j * b_dim1] = b[i__ + j * b_dim1] - dl[
				i__] * temp;
#line 265 "sgtts2.f"
			b[i__ + j * b_dim1] = temp;
#line 266 "sgtts2.f"
		    }
#line 267 "sgtts2.f"
/* L110: */
#line 267 "sgtts2.f"
		}
#line 268 "sgtts2.f"
/* L120: */
#line 268 "sgtts2.f"
	    }
#line 269 "sgtts2.f"
	}
#line 270 "sgtts2.f"
    }

/*     End of SGTTS2 */

#line 274 "sgtts2.f"
    return 0;
} /* sgtts2_ */

