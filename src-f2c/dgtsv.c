#line 1 "dgtsv.f"
/* dgtsv.f -- translated by f2c (version 20100827).
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

#line 1 "dgtsv.f"
/* > \brief <b> DGTSV computes the solution to system of linear equations A * X = B for GT matrices <b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGTSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtsv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtsv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtsv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGTSV  solves the equation */
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
/* >          DL is DOUBLE PRECISION array, dimension (N-1) */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, D must contain the diagonal elements of A. */
/* > */
/* >          On exit, D is overwritten by the n diagonal elements of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* >          DU is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, DU must contain the (n-1) super-diagonal elements */
/* >          of A. */
/* > */
/* >          On exit, DU is overwritten by the (n-1) elements of the first */
/* >          super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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

/* > \ingroup doubleGTsolve */

/*  ===================================================================== */
/* Subroutine */ int dgtsv_(integer *n, integer *nrhs, doublereal *dl, 
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

#line 160 "dgtsv.f"
    /* Parameter adjustments */
#line 160 "dgtsv.f"
    --dl;
#line 160 "dgtsv.f"
    --d__;
#line 160 "dgtsv.f"
    --du;
#line 160 "dgtsv.f"
    b_dim1 = *ldb;
#line 160 "dgtsv.f"
    b_offset = 1 + b_dim1;
#line 160 "dgtsv.f"
    b -= b_offset;
#line 160 "dgtsv.f"

#line 160 "dgtsv.f"
    /* Function Body */
#line 160 "dgtsv.f"
    *info = 0;
#line 161 "dgtsv.f"
    if (*n < 0) {
#line 162 "dgtsv.f"
	*info = -1;
#line 163 "dgtsv.f"
    } else if (*nrhs < 0) {
#line 164 "dgtsv.f"
	*info = -2;
#line 165 "dgtsv.f"
    } else if (*ldb < max(1,*n)) {
#line 166 "dgtsv.f"
	*info = -7;
#line 167 "dgtsv.f"
    }
#line 168 "dgtsv.f"
    if (*info != 0) {
#line 169 "dgtsv.f"
	i__1 = -(*info);
#line 169 "dgtsv.f"
	xerbla_("DGTSV ", &i__1, (ftnlen)6);
#line 170 "dgtsv.f"
	return 0;
#line 171 "dgtsv.f"
    }

#line 173 "dgtsv.f"
    if (*n == 0) {
#line 173 "dgtsv.f"
	return 0;
#line 173 "dgtsv.f"
    }

#line 176 "dgtsv.f"
    if (*nrhs == 1) {
#line 177 "dgtsv.f"
	i__1 = *n - 2;
#line 177 "dgtsv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 178 "dgtsv.f"
	    if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {

/*              No row interchange required */

#line 182 "dgtsv.f"
		if (d__[i__] != 0.) {
#line 183 "dgtsv.f"
		    fact = dl[i__] / d__[i__];
#line 184 "dgtsv.f"
		    d__[i__ + 1] -= fact * du[i__];
#line 185 "dgtsv.f"
		    b[i__ + 1 + b_dim1] -= fact * b[i__ + b_dim1];
#line 186 "dgtsv.f"
		} else {
#line 187 "dgtsv.f"
		    *info = i__;
#line 188 "dgtsv.f"
		    return 0;
#line 189 "dgtsv.f"
		}
#line 190 "dgtsv.f"
		dl[i__] = 0.;
#line 191 "dgtsv.f"
	    } else {

/*              Interchange rows I and I+1 */

#line 195 "dgtsv.f"
		fact = d__[i__] / dl[i__];
#line 196 "dgtsv.f"
		d__[i__] = dl[i__];
#line 197 "dgtsv.f"
		temp = d__[i__ + 1];
#line 198 "dgtsv.f"
		d__[i__ + 1] = du[i__] - fact * temp;
#line 199 "dgtsv.f"
		dl[i__] = du[i__ + 1];
#line 200 "dgtsv.f"
		du[i__ + 1] = -fact * dl[i__];
#line 201 "dgtsv.f"
		du[i__] = temp;
#line 202 "dgtsv.f"
		temp = b[i__ + b_dim1];
#line 203 "dgtsv.f"
		b[i__ + b_dim1] = b[i__ + 1 + b_dim1];
#line 204 "dgtsv.f"
		b[i__ + 1 + b_dim1] = temp - fact * b[i__ + 1 + b_dim1];
#line 205 "dgtsv.f"
	    }
#line 206 "dgtsv.f"
/* L10: */
#line 206 "dgtsv.f"
	}
#line 207 "dgtsv.f"
	if (*n > 1) {
#line 208 "dgtsv.f"
	    i__ = *n - 1;
#line 209 "dgtsv.f"
	    if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {
#line 210 "dgtsv.f"
		if (d__[i__] != 0.) {
#line 211 "dgtsv.f"
		    fact = dl[i__] / d__[i__];
#line 212 "dgtsv.f"
		    d__[i__ + 1] -= fact * du[i__];
#line 213 "dgtsv.f"
		    b[i__ + 1 + b_dim1] -= fact * b[i__ + b_dim1];
#line 214 "dgtsv.f"
		} else {
#line 215 "dgtsv.f"
		    *info = i__;
#line 216 "dgtsv.f"
		    return 0;
#line 217 "dgtsv.f"
		}
#line 218 "dgtsv.f"
	    } else {
#line 219 "dgtsv.f"
		fact = d__[i__] / dl[i__];
#line 220 "dgtsv.f"
		d__[i__] = dl[i__];
#line 221 "dgtsv.f"
		temp = d__[i__ + 1];
#line 222 "dgtsv.f"
		d__[i__ + 1] = du[i__] - fact * temp;
#line 223 "dgtsv.f"
		du[i__] = temp;
#line 224 "dgtsv.f"
		temp = b[i__ + b_dim1];
#line 225 "dgtsv.f"
		b[i__ + b_dim1] = b[i__ + 1 + b_dim1];
#line 226 "dgtsv.f"
		b[i__ + 1 + b_dim1] = temp - fact * b[i__ + 1 + b_dim1];
#line 227 "dgtsv.f"
	    }
#line 228 "dgtsv.f"
	}
#line 229 "dgtsv.f"
	if (d__[*n] == 0.) {
#line 230 "dgtsv.f"
	    *info = *n;
#line 231 "dgtsv.f"
	    return 0;
#line 232 "dgtsv.f"
	}
#line 233 "dgtsv.f"
    } else {
#line 234 "dgtsv.f"
	i__1 = *n - 2;
#line 234 "dgtsv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "dgtsv.f"
	    if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {

/*              No row interchange required */

#line 239 "dgtsv.f"
		if (d__[i__] != 0.) {
#line 240 "dgtsv.f"
		    fact = dl[i__] / d__[i__];
#line 241 "dgtsv.f"
		    d__[i__ + 1] -= fact * du[i__];
#line 242 "dgtsv.f"
		    i__2 = *nrhs;
#line 242 "dgtsv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 243 "dgtsv.f"
			b[i__ + 1 + j * b_dim1] -= fact * b[i__ + j * b_dim1];
#line 244 "dgtsv.f"
/* L20: */
#line 244 "dgtsv.f"
		    }
#line 245 "dgtsv.f"
		} else {
#line 246 "dgtsv.f"
		    *info = i__;
#line 247 "dgtsv.f"
		    return 0;
#line 248 "dgtsv.f"
		}
#line 249 "dgtsv.f"
		dl[i__] = 0.;
#line 250 "dgtsv.f"
	    } else {

/*              Interchange rows I and I+1 */

#line 254 "dgtsv.f"
		fact = d__[i__] / dl[i__];
#line 255 "dgtsv.f"
		d__[i__] = dl[i__];
#line 256 "dgtsv.f"
		temp = d__[i__ + 1];
#line 257 "dgtsv.f"
		d__[i__ + 1] = du[i__] - fact * temp;
#line 258 "dgtsv.f"
		dl[i__] = du[i__ + 1];
#line 259 "dgtsv.f"
		du[i__ + 1] = -fact * dl[i__];
#line 260 "dgtsv.f"
		du[i__] = temp;
#line 261 "dgtsv.f"
		i__2 = *nrhs;
#line 261 "dgtsv.f"
		for (j = 1; j <= i__2; ++j) {
#line 262 "dgtsv.f"
		    temp = b[i__ + j * b_dim1];
#line 263 "dgtsv.f"
		    b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
#line 264 "dgtsv.f"
		    b[i__ + 1 + j * b_dim1] = temp - fact * b[i__ + 1 + j * 
			    b_dim1];
#line 265 "dgtsv.f"
/* L30: */
#line 265 "dgtsv.f"
		}
#line 266 "dgtsv.f"
	    }
#line 267 "dgtsv.f"
/* L40: */
#line 267 "dgtsv.f"
	}
#line 268 "dgtsv.f"
	if (*n > 1) {
#line 269 "dgtsv.f"
	    i__ = *n - 1;
#line 270 "dgtsv.f"
	    if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {
#line 271 "dgtsv.f"
		if (d__[i__] != 0.) {
#line 272 "dgtsv.f"
		    fact = dl[i__] / d__[i__];
#line 273 "dgtsv.f"
		    d__[i__ + 1] -= fact * du[i__];
#line 274 "dgtsv.f"
		    i__1 = *nrhs;
#line 274 "dgtsv.f"
		    for (j = 1; j <= i__1; ++j) {
#line 275 "dgtsv.f"
			b[i__ + 1 + j * b_dim1] -= fact * b[i__ + j * b_dim1];
#line 276 "dgtsv.f"
/* L50: */
#line 276 "dgtsv.f"
		    }
#line 277 "dgtsv.f"
		} else {
#line 278 "dgtsv.f"
		    *info = i__;
#line 279 "dgtsv.f"
		    return 0;
#line 280 "dgtsv.f"
		}
#line 281 "dgtsv.f"
	    } else {
#line 282 "dgtsv.f"
		fact = d__[i__] / dl[i__];
#line 283 "dgtsv.f"
		d__[i__] = dl[i__];
#line 284 "dgtsv.f"
		temp = d__[i__ + 1];
#line 285 "dgtsv.f"
		d__[i__ + 1] = du[i__] - fact * temp;
#line 286 "dgtsv.f"
		du[i__] = temp;
#line 287 "dgtsv.f"
		i__1 = *nrhs;
#line 287 "dgtsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 288 "dgtsv.f"
		    temp = b[i__ + j * b_dim1];
#line 289 "dgtsv.f"
		    b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
#line 290 "dgtsv.f"
		    b[i__ + 1 + j * b_dim1] = temp - fact * b[i__ + 1 + j * 
			    b_dim1];
#line 291 "dgtsv.f"
/* L60: */
#line 291 "dgtsv.f"
		}
#line 292 "dgtsv.f"
	    }
#line 293 "dgtsv.f"
	}
#line 294 "dgtsv.f"
	if (d__[*n] == 0.) {
#line 295 "dgtsv.f"
	    *info = *n;
#line 296 "dgtsv.f"
	    return 0;
#line 297 "dgtsv.f"
	}
#line 298 "dgtsv.f"
    }

/*     Back solve with the matrix U from the factorization. */

#line 302 "dgtsv.f"
    if (*nrhs <= 2) {
#line 303 "dgtsv.f"
	j = 1;
#line 304 "dgtsv.f"
L70:
#line 305 "dgtsv.f"
	b[*n + j * b_dim1] /= d__[*n];
#line 306 "dgtsv.f"
	if (*n > 1) {
#line 306 "dgtsv.f"
	    b[*n - 1 + j * b_dim1] = (b[*n - 1 + j * b_dim1] - du[*n - 1] * b[
		    *n + j * b_dim1]) / d__[*n - 1];
#line 306 "dgtsv.f"
	}
#line 308 "dgtsv.f"
	for (i__ = *n - 2; i__ >= 1; --i__) {
#line 309 "dgtsv.f"
	    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ + 1 
		    + j * b_dim1] - dl[i__] * b[i__ + 2 + j * b_dim1]) / d__[
		    i__];
#line 311 "dgtsv.f"
/* L80: */
#line 311 "dgtsv.f"
	}
#line 312 "dgtsv.f"
	if (j < *nrhs) {
#line 313 "dgtsv.f"
	    ++j;
#line 314 "dgtsv.f"
	    goto L70;
#line 315 "dgtsv.f"
	}
#line 316 "dgtsv.f"
    } else {
#line 317 "dgtsv.f"
	i__1 = *nrhs;
#line 317 "dgtsv.f"
	for (j = 1; j <= i__1; ++j) {
#line 318 "dgtsv.f"
	    b[*n + j * b_dim1] /= d__[*n];
#line 319 "dgtsv.f"
	    if (*n > 1) {
#line 319 "dgtsv.f"
		b[*n - 1 + j * b_dim1] = (b[*n - 1 + j * b_dim1] - du[*n - 1] 
			* b[*n + j * b_dim1]) / d__[*n - 1];
#line 319 "dgtsv.f"
	    }
#line 322 "dgtsv.f"
	    for (i__ = *n - 2; i__ >= 1; --i__) {
#line 323 "dgtsv.f"
		b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ 
			+ 1 + j * b_dim1] - dl[i__] * b[i__ + 2 + j * b_dim1])
			 / d__[i__];
#line 325 "dgtsv.f"
/* L90: */
#line 325 "dgtsv.f"
	    }
#line 326 "dgtsv.f"
/* L100: */
#line 326 "dgtsv.f"
	}
#line 327 "dgtsv.f"
    }

#line 329 "dgtsv.f"
    return 0;

/*     End of DGTSV */

} /* dgtsv_ */

