#line 1 "sgtsv.f"
/* sgtsv.f -- translated by f2c (version 20100827).
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

#line 1 "sgtsv.f"
/* > \brief <b> SGTSV computes the solution to system of linear equations A * X = B for GT matrices <b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGTSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtsv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtsv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtsv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGTSV( N, NRHS, DL, D, DU, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               B( LDB, * ), D( * ), DL( * ), DU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGTSV  solves the equation */
/* > */
/* >    A*X = B, */
/* > */
/* > where A is an n by n tridiagonal matrix, by Gaussian elimination with */
/* > partial pivoting. */
/* > */
/* > Note that the equation  A**T*X = B  may be solved by interchanging the */
/* > order of the arguments DU and DL. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* >          DL is REAL array, dimension (N-1) */
/* >          On entry, DL must contain the (n-1) sub-diagonal elements of */
/* >          A. */
/* > */
/* >          On exit, DL is overwritten by the (n-2) elements of the */
/* >          second super-diagonal of the upper triangular matrix U from */
/* >          the LU factorization of A, in DL(1), ..., DL(n-2). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, D must contain the diagonal elements of A. */
/* > */
/* >          On exit, D is overwritten by the n diagonal elements of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* >          DU is REAL array, dimension (N-1) */
/* >          On entry, DU must contain the (n-1) super-diagonal elements */
/* >          of A. */
/* > */
/* >          On exit, DU is overwritten by the (n-1) elements of the first */
/* >          super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
/* >          On entry, the N by NRHS matrix of right hand side matrix B. */
/* >          On exit, if INFO = 0, the N by NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, U(i,i) is exactly zero, and the solution */
/* >               has not been computed.  The factorization has not been */
/* >               completed unless i = N. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGTsolve */

/*  ===================================================================== */
/* Subroutine */ int sgtsv_(integer *n, integer *nrhs, doublereal *dl, 
	doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer 
	*info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal fact, temp;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.2) -- */
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

#line 160 "sgtsv.f"
    /* Parameter adjustments */
#line 160 "sgtsv.f"
    --dl;
#line 160 "sgtsv.f"
    --d__;
#line 160 "sgtsv.f"
    --du;
#line 160 "sgtsv.f"
    b_dim1 = *ldb;
#line 160 "sgtsv.f"
    b_offset = 1 + b_dim1;
#line 160 "sgtsv.f"
    b -= b_offset;
#line 160 "sgtsv.f"

#line 160 "sgtsv.f"
    /* Function Body */
#line 160 "sgtsv.f"
    *info = 0;
#line 161 "sgtsv.f"
    if (*n < 0) {
#line 162 "sgtsv.f"
	*info = -1;
#line 163 "sgtsv.f"
    } else if (*nrhs < 0) {
#line 164 "sgtsv.f"
	*info = -2;
#line 165 "sgtsv.f"
    } else if (*ldb < max(1,*n)) {
#line 166 "sgtsv.f"
	*info = -7;
#line 167 "sgtsv.f"
    }
#line 168 "sgtsv.f"
    if (*info != 0) {
#line 169 "sgtsv.f"
	i__1 = -(*info);
#line 169 "sgtsv.f"
	xerbla_("SGTSV ", &i__1, (ftnlen)6);
#line 170 "sgtsv.f"
	return 0;
#line 171 "sgtsv.f"
    }

#line 173 "sgtsv.f"
    if (*n == 0) {
#line 173 "sgtsv.f"
	return 0;
#line 173 "sgtsv.f"
    }

#line 176 "sgtsv.f"
    if (*nrhs == 1) {
#line 177 "sgtsv.f"
	i__1 = *n - 2;
#line 177 "sgtsv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 178 "sgtsv.f"
	    if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {

/*              No row interchange required */

#line 182 "sgtsv.f"
		if (d__[i__] != 0.) {
#line 183 "sgtsv.f"
		    fact = dl[i__] / d__[i__];
#line 184 "sgtsv.f"
		    d__[i__ + 1] -= fact * du[i__];
#line 185 "sgtsv.f"
		    b[i__ + 1 + b_dim1] -= fact * b[i__ + b_dim1];
#line 186 "sgtsv.f"
		} else {
#line 187 "sgtsv.f"
		    *info = i__;
#line 188 "sgtsv.f"
		    return 0;
#line 189 "sgtsv.f"
		}
#line 190 "sgtsv.f"
		dl[i__] = 0.;
#line 191 "sgtsv.f"
	    } else {

/*              Interchange rows I and I+1 */

#line 195 "sgtsv.f"
		fact = d__[i__] / dl[i__];
#line 196 "sgtsv.f"
		d__[i__] = dl[i__];
#line 197 "sgtsv.f"
		temp = d__[i__ + 1];
#line 198 "sgtsv.f"
		d__[i__ + 1] = du[i__] - fact * temp;
#line 199 "sgtsv.f"
		dl[i__] = du[i__ + 1];
#line 200 "sgtsv.f"
		du[i__ + 1] = -fact * dl[i__];
#line 201 "sgtsv.f"
		du[i__] = temp;
#line 202 "sgtsv.f"
		temp = b[i__ + b_dim1];
#line 203 "sgtsv.f"
		b[i__ + b_dim1] = b[i__ + 1 + b_dim1];
#line 204 "sgtsv.f"
		b[i__ + 1 + b_dim1] = temp - fact * b[i__ + 1 + b_dim1];
#line 205 "sgtsv.f"
	    }
#line 206 "sgtsv.f"
/* L10: */
#line 206 "sgtsv.f"
	}
#line 207 "sgtsv.f"
	if (*n > 1) {
#line 208 "sgtsv.f"
	    i__ = *n - 1;
#line 209 "sgtsv.f"
	    if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {
#line 210 "sgtsv.f"
		if (d__[i__] != 0.) {
#line 211 "sgtsv.f"
		    fact = dl[i__] / d__[i__];
#line 212 "sgtsv.f"
		    d__[i__ + 1] -= fact * du[i__];
#line 213 "sgtsv.f"
		    b[i__ + 1 + b_dim1] -= fact * b[i__ + b_dim1];
#line 214 "sgtsv.f"
		} else {
#line 215 "sgtsv.f"
		    *info = i__;
#line 216 "sgtsv.f"
		    return 0;
#line 217 "sgtsv.f"
		}
#line 218 "sgtsv.f"
	    } else {
#line 219 "sgtsv.f"
		fact = d__[i__] / dl[i__];
#line 220 "sgtsv.f"
		d__[i__] = dl[i__];
#line 221 "sgtsv.f"
		temp = d__[i__ + 1];
#line 222 "sgtsv.f"
		d__[i__ + 1] = du[i__] - fact * temp;
#line 223 "sgtsv.f"
		du[i__] = temp;
#line 224 "sgtsv.f"
		temp = b[i__ + b_dim1];
#line 225 "sgtsv.f"
		b[i__ + b_dim1] = b[i__ + 1 + b_dim1];
#line 226 "sgtsv.f"
		b[i__ + 1 + b_dim1] = temp - fact * b[i__ + 1 + b_dim1];
#line 227 "sgtsv.f"
	    }
#line 228 "sgtsv.f"
	}
#line 229 "sgtsv.f"
	if (d__[*n] == 0.) {
#line 230 "sgtsv.f"
	    *info = *n;
#line 231 "sgtsv.f"
	    return 0;
#line 232 "sgtsv.f"
	}
#line 233 "sgtsv.f"
    } else {
#line 234 "sgtsv.f"
	i__1 = *n - 2;
#line 234 "sgtsv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "sgtsv.f"
	    if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {

/*              No row interchange required */

#line 239 "sgtsv.f"
		if (d__[i__] != 0.) {
#line 240 "sgtsv.f"
		    fact = dl[i__] / d__[i__];
#line 241 "sgtsv.f"
		    d__[i__ + 1] -= fact * du[i__];
#line 242 "sgtsv.f"
		    i__2 = *nrhs;
#line 242 "sgtsv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 243 "sgtsv.f"
			b[i__ + 1 + j * b_dim1] -= fact * b[i__ + j * b_dim1];
#line 244 "sgtsv.f"
/* L20: */
#line 244 "sgtsv.f"
		    }
#line 245 "sgtsv.f"
		} else {
#line 246 "sgtsv.f"
		    *info = i__;
#line 247 "sgtsv.f"
		    return 0;
#line 248 "sgtsv.f"
		}
#line 249 "sgtsv.f"
		dl[i__] = 0.;
#line 250 "sgtsv.f"
	    } else {

/*              Interchange rows I and I+1 */

#line 254 "sgtsv.f"
		fact = d__[i__] / dl[i__];
#line 255 "sgtsv.f"
		d__[i__] = dl[i__];
#line 256 "sgtsv.f"
		temp = d__[i__ + 1];
#line 257 "sgtsv.f"
		d__[i__ + 1] = du[i__] - fact * temp;
#line 258 "sgtsv.f"
		dl[i__] = du[i__ + 1];
#line 259 "sgtsv.f"
		du[i__ + 1] = -fact * dl[i__];
#line 260 "sgtsv.f"
		du[i__] = temp;
#line 261 "sgtsv.f"
		i__2 = *nrhs;
#line 261 "sgtsv.f"
		for (j = 1; j <= i__2; ++j) {
#line 262 "sgtsv.f"
		    temp = b[i__ + j * b_dim1];
#line 263 "sgtsv.f"
		    b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
#line 264 "sgtsv.f"
		    b[i__ + 1 + j * b_dim1] = temp - fact * b[i__ + 1 + j * 
			    b_dim1];
#line 265 "sgtsv.f"
/* L30: */
#line 265 "sgtsv.f"
		}
#line 266 "sgtsv.f"
	    }
#line 267 "sgtsv.f"
/* L40: */
#line 267 "sgtsv.f"
	}
#line 268 "sgtsv.f"
	if (*n > 1) {
#line 269 "sgtsv.f"
	    i__ = *n - 1;
#line 270 "sgtsv.f"
	    if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {
#line 271 "sgtsv.f"
		if (d__[i__] != 0.) {
#line 272 "sgtsv.f"
		    fact = dl[i__] / d__[i__];
#line 273 "sgtsv.f"
		    d__[i__ + 1] -= fact * du[i__];
#line 274 "sgtsv.f"
		    i__1 = *nrhs;
#line 274 "sgtsv.f"
		    for (j = 1; j <= i__1; ++j) {
#line 275 "sgtsv.f"
			b[i__ + 1 + j * b_dim1] -= fact * b[i__ + j * b_dim1];
#line 276 "sgtsv.f"
/* L50: */
#line 276 "sgtsv.f"
		    }
#line 277 "sgtsv.f"
		} else {
#line 278 "sgtsv.f"
		    *info = i__;
#line 279 "sgtsv.f"
		    return 0;
#line 280 "sgtsv.f"
		}
#line 281 "sgtsv.f"
	    } else {
#line 282 "sgtsv.f"
		fact = d__[i__] / dl[i__];
#line 283 "sgtsv.f"
		d__[i__] = dl[i__];
#line 284 "sgtsv.f"
		temp = d__[i__ + 1];
#line 285 "sgtsv.f"
		d__[i__ + 1] = du[i__] - fact * temp;
#line 286 "sgtsv.f"
		du[i__] = temp;
#line 287 "sgtsv.f"
		i__1 = *nrhs;
#line 287 "sgtsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 288 "sgtsv.f"
		    temp = b[i__ + j * b_dim1];
#line 289 "sgtsv.f"
		    b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
#line 290 "sgtsv.f"
		    b[i__ + 1 + j * b_dim1] = temp - fact * b[i__ + 1 + j * 
			    b_dim1];
#line 291 "sgtsv.f"
/* L60: */
#line 291 "sgtsv.f"
		}
#line 292 "sgtsv.f"
	    }
#line 293 "sgtsv.f"
	}
#line 294 "sgtsv.f"
	if (d__[*n] == 0.) {
#line 295 "sgtsv.f"
	    *info = *n;
#line 296 "sgtsv.f"
	    return 0;
#line 297 "sgtsv.f"
	}
#line 298 "sgtsv.f"
    }

/*     Back solve with the matrix U from the factorization. */

#line 302 "sgtsv.f"
    if (*nrhs <= 2) {
#line 303 "sgtsv.f"
	j = 1;
#line 304 "sgtsv.f"
L70:
#line 305 "sgtsv.f"
	b[*n + j * b_dim1] /= d__[*n];
#line 306 "sgtsv.f"
	if (*n > 1) {
#line 306 "sgtsv.f"
	    b[*n - 1 + j * b_dim1] = (b[*n - 1 + j * b_dim1] - du[*n - 1] * b[
		    *n + j * b_dim1]) / d__[*n - 1];
#line 306 "sgtsv.f"
	}
#line 308 "sgtsv.f"
	for (i__ = *n - 2; i__ >= 1; --i__) {
#line 309 "sgtsv.f"
	    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ + 1 
		    + j * b_dim1] - dl[i__] * b[i__ + 2 + j * b_dim1]) / d__[
		    i__];
#line 311 "sgtsv.f"
/* L80: */
#line 311 "sgtsv.f"
	}
#line 312 "sgtsv.f"
	if (j < *nrhs) {
#line 313 "sgtsv.f"
	    ++j;
#line 314 "sgtsv.f"
	    goto L70;
#line 315 "sgtsv.f"
	}
#line 316 "sgtsv.f"
    } else {
#line 317 "sgtsv.f"
	i__1 = *nrhs;
#line 317 "sgtsv.f"
	for (j = 1; j <= i__1; ++j) {
#line 318 "sgtsv.f"
	    b[*n + j * b_dim1] /= d__[*n];
#line 319 "sgtsv.f"
	    if (*n > 1) {
#line 319 "sgtsv.f"
		b[*n - 1 + j * b_dim1] = (b[*n - 1 + j * b_dim1] - du[*n - 1] 
			* b[*n + j * b_dim1]) / d__[*n - 1];
#line 319 "sgtsv.f"
	    }
#line 322 "sgtsv.f"
	    for (i__ = *n - 2; i__ >= 1; --i__) {
#line 323 "sgtsv.f"
		b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ 
			+ 1 + j * b_dim1] - dl[i__] * b[i__ + 2 + j * b_dim1])
			 / d__[i__];
#line 325 "sgtsv.f"
/* L90: */
#line 325 "sgtsv.f"
	    }
#line 326 "sgtsv.f"
/* L100: */
#line 326 "sgtsv.f"
	}
#line 327 "sgtsv.f"
    }

#line 329 "sgtsv.f"
    return 0;

/*     End of SGTSV */

} /* sgtsv_ */

