#line 1 "dlagts.f"
/* dlagts.f -- translated by f2c (version 20100827).
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

#line 1 "dlagts.f"
/* > \brief \b DLAGTS solves the system of equations (T-λI)x = y or (T-λI)Tx = y,where T is a general tridia
gonal matrix and λ a scalar, using the LU factorization computed by slagtf. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAGTS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagts.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagts.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagts.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, JOB, N */
/*       DOUBLE PRECISION   TOL */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IN( * ) */
/*       DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAGTS may be used to solve one of the systems of equations */
/* > */
/* >    (T - lambda*I)*x = y   or   (T - lambda*I)**T*x = y, */
/* > */
/* > where T is an n by n tridiagonal matrix, for x, following the */
/* > factorization of (T - lambda*I) as */
/* > */
/* >    (T - lambda*I) = P*L*U , */
/* > */
/* > by routine DLAGTF. The choice of equation to be solved is */
/* > controlled by the argument JOB, and in each case there is an option */
/* > to perturb zero or very small diagonal elements of U, this option */
/* > being intended for use in applications such as inverse iteration. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is INTEGER */
/* >          Specifies the job to be performed by DLAGTS as follows: */
/* >          =  1: The equations  (T - lambda*I)x = y  are to be solved, */
/* >                but diagonal elements of U are not to be perturbed. */
/* >          = -1: The equations  (T - lambda*I)x = y  are to be solved */
/* >                and, if overflow would otherwise occur, the diagonal */
/* >                elements of U are to be perturbed. See argument TOL */
/* >                below. */
/* >          =  2: The equations  (T - lambda*I)**Tx = y  are to be solved, */
/* >                but diagonal elements of U are not to be perturbed. */
/* >          = -2: The equations  (T - lambda*I)**Tx = y  are to be solved */
/* >                and, if overflow would otherwise occur, the diagonal */
/* >                elements of U are to be perturbed. See argument TOL */
/* >                below. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, A must contain the diagonal elements of U as */
/* >          returned from DLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, B must contain the first super-diagonal elements of */
/* >          U as returned from DLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, C must contain the sub-diagonal elements of L as */
/* >          returned from DLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N-2) */
/* >          On entry, D must contain the second super-diagonal elements */
/* >          of U as returned from DLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in] IN */
/* > \verbatim */
/* >          IN is INTEGER array, dimension (N) */
/* >          On entry, IN must contain details of the matrix P as returned */
/* >          from DLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the right hand side vector y. */
/* >          On exit, Y is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[in,out] TOL */
/* > \verbatim */
/* >          TOL is DOUBLE PRECISION */
/* >          On entry, with  JOB .lt. 0, TOL should be the minimum */
/* >          perturbation to be made to very small diagonal elements of U. */
/* >          TOL should normally be chosen as about eps*norm(U), where eps */
/* >          is the relative machine precision, but if TOL is supplied as */
/* >          non-positive, then it is reset to eps*max( abs( u(i,j) ) ). */
/* >          If  JOB .gt. 0  then TOL is not referenced. */
/* > */
/* >          On exit, TOL is changed as described above, only if TOL is */
/* >          non-positive on entry. Otherwise TOL is unchanged. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0   : successful exit */
/* >          .lt. 0: if INFO = -i, the i-th argument had an illegal value */
/* >          .gt. 0: overflow would occur when computing the INFO(th) */
/* >                  element of the solution vector x. This can only occur */
/* >                  when JOB is supplied as positive and either means */
/* >                  that a diagonal element of U is very small, or that */
/* >                  the elements of the right-hand side vector y are very */
/* >                  large. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlagts_(integer *job, integer *n, doublereal *a, 
	doublereal *b, doublereal *c__, doublereal *d__, integer *in, 
	doublereal *y, doublereal *tol, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer k;
    static doublereal ak, eps, temp, pert, absak, sfmin;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

#line 200 "dlagts.f"
    /* Parameter adjustments */
#line 200 "dlagts.f"
    --y;
#line 200 "dlagts.f"
    --in;
#line 200 "dlagts.f"
    --d__;
#line 200 "dlagts.f"
    --c__;
#line 200 "dlagts.f"
    --b;
#line 200 "dlagts.f"
    --a;
#line 200 "dlagts.f"

#line 200 "dlagts.f"
    /* Function Body */
#line 200 "dlagts.f"
    *info = 0;
#line 201 "dlagts.f"
    if (abs(*job) > 2 || *job == 0) {
#line 202 "dlagts.f"
	*info = -1;
#line 203 "dlagts.f"
    } else if (*n < 0) {
#line 204 "dlagts.f"
	*info = -2;
#line 205 "dlagts.f"
    }
#line 206 "dlagts.f"
    if (*info != 0) {
#line 207 "dlagts.f"
	i__1 = -(*info);
#line 207 "dlagts.f"
	xerbla_("DLAGTS", &i__1, (ftnlen)6);
#line 208 "dlagts.f"
	return 0;
#line 209 "dlagts.f"
    }

#line 211 "dlagts.f"
    if (*n == 0) {
#line 211 "dlagts.f"
	return 0;
#line 211 "dlagts.f"
    }

#line 214 "dlagts.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 215 "dlagts.f"
    sfmin = dlamch_("Safe minimum", (ftnlen)12);
#line 216 "dlagts.f"
    bignum = 1. / sfmin;

#line 218 "dlagts.f"
    if (*job < 0) {
#line 219 "dlagts.f"
	if (*tol <= 0.) {
#line 220 "dlagts.f"
	    *tol = abs(a[1]);
#line 221 "dlagts.f"
	    if (*n > 1) {
/* Computing MAX */
#line 221 "dlagts.f"
		d__1 = *tol, d__2 = abs(a[2]), d__1 = max(d__1,d__2), d__2 = 
			abs(b[1]);
#line 221 "dlagts.f"
		*tol = max(d__1,d__2);
#line 221 "dlagts.f"
	    }
#line 223 "dlagts.f"
	    i__1 = *n;
#line 223 "dlagts.f"
	    for (k = 3; k <= i__1; ++k) {
/* Computing MAX */
#line 224 "dlagts.f"
		d__4 = *tol, d__5 = (d__1 = a[k], abs(d__1)), d__4 = max(d__4,
			d__5), d__5 = (d__2 = b[k - 1], abs(d__2)), d__4 = 
			max(d__4,d__5), d__5 = (d__3 = d__[k - 2], abs(d__3));
#line 224 "dlagts.f"
		*tol = max(d__4,d__5);
#line 226 "dlagts.f"
/* L10: */
#line 226 "dlagts.f"
	    }
#line 227 "dlagts.f"
	    *tol *= eps;
#line 228 "dlagts.f"
	    if (*tol == 0.) {
#line 228 "dlagts.f"
		*tol = eps;
#line 228 "dlagts.f"
	    }
#line 230 "dlagts.f"
	}
#line 231 "dlagts.f"
    }

#line 233 "dlagts.f"
    if (abs(*job) == 1) {
#line 234 "dlagts.f"
	i__1 = *n;
#line 234 "dlagts.f"
	for (k = 2; k <= i__1; ++k) {
#line 235 "dlagts.f"
	    if (in[k - 1] == 0) {
#line 236 "dlagts.f"
		y[k] -= c__[k - 1] * y[k - 1];
#line 237 "dlagts.f"
	    } else {
#line 238 "dlagts.f"
		temp = y[k - 1];
#line 239 "dlagts.f"
		y[k - 1] = y[k];
#line 240 "dlagts.f"
		y[k] = temp - c__[k - 1] * y[k];
#line 241 "dlagts.f"
	    }
#line 242 "dlagts.f"
/* L20: */
#line 242 "dlagts.f"
	}
#line 243 "dlagts.f"
	if (*job == 1) {
#line 244 "dlagts.f"
	    for (k = *n; k >= 1; --k) {
#line 245 "dlagts.f"
		if (k <= *n - 2) {
#line 246 "dlagts.f"
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
#line 247 "dlagts.f"
		} else if (k == *n - 1) {
#line 248 "dlagts.f"
		    temp = y[k] - b[k] * y[k + 1];
#line 249 "dlagts.f"
		} else {
#line 250 "dlagts.f"
		    temp = y[k];
#line 251 "dlagts.f"
		}
#line 252 "dlagts.f"
		ak = a[k];
#line 253 "dlagts.f"
		absak = abs(ak);
#line 254 "dlagts.f"
		if (absak < 1.) {
#line 255 "dlagts.f"
		    if (absak < sfmin) {
#line 256 "dlagts.f"
			if (absak == 0. || abs(temp) * sfmin > absak) {
#line 258 "dlagts.f"
			    *info = k;
#line 259 "dlagts.f"
			    return 0;
#line 260 "dlagts.f"
			} else {
#line 261 "dlagts.f"
			    temp *= bignum;
#line 262 "dlagts.f"
			    ak *= bignum;
#line 263 "dlagts.f"
			}
#line 264 "dlagts.f"
		    } else if (abs(temp) > absak * bignum) {
#line 265 "dlagts.f"
			*info = k;
#line 266 "dlagts.f"
			return 0;
#line 267 "dlagts.f"
		    }
#line 268 "dlagts.f"
		}
#line 269 "dlagts.f"
		y[k] = temp / ak;
#line 270 "dlagts.f"
/* L30: */
#line 270 "dlagts.f"
	    }
#line 271 "dlagts.f"
	} else {
#line 272 "dlagts.f"
	    for (k = *n; k >= 1; --k) {
#line 273 "dlagts.f"
		if (k <= *n - 2) {
#line 274 "dlagts.f"
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
#line 275 "dlagts.f"
		} else if (k == *n - 1) {
#line 276 "dlagts.f"
		    temp = y[k] - b[k] * y[k + 1];
#line 277 "dlagts.f"
		} else {
#line 278 "dlagts.f"
		    temp = y[k];
#line 279 "dlagts.f"
		}
#line 280 "dlagts.f"
		ak = a[k];
#line 281 "dlagts.f"
		pert = d_sign(tol, &ak);
#line 282 "dlagts.f"
L40:
#line 283 "dlagts.f"
		absak = abs(ak);
#line 284 "dlagts.f"
		if (absak < 1.) {
#line 285 "dlagts.f"
		    if (absak < sfmin) {
#line 286 "dlagts.f"
			if (absak == 0. || abs(temp) * sfmin > absak) {
#line 288 "dlagts.f"
			    ak += pert;
#line 289 "dlagts.f"
			    pert *= 2;
#line 290 "dlagts.f"
			    goto L40;
#line 291 "dlagts.f"
			} else {
#line 292 "dlagts.f"
			    temp *= bignum;
#line 293 "dlagts.f"
			    ak *= bignum;
#line 294 "dlagts.f"
			}
#line 295 "dlagts.f"
		    } else if (abs(temp) > absak * bignum) {
#line 296 "dlagts.f"
			ak += pert;
#line 297 "dlagts.f"
			pert *= 2;
#line 298 "dlagts.f"
			goto L40;
#line 299 "dlagts.f"
		    }
#line 300 "dlagts.f"
		}
#line 301 "dlagts.f"
		y[k] = temp / ak;
#line 302 "dlagts.f"
/* L50: */
#line 302 "dlagts.f"
	    }
#line 303 "dlagts.f"
	}
#line 304 "dlagts.f"
    } else {

/*        Come to here if  JOB = 2 or -2 */

#line 308 "dlagts.f"
	if (*job == 2) {
#line 309 "dlagts.f"
	    i__1 = *n;
#line 309 "dlagts.f"
	    for (k = 1; k <= i__1; ++k) {
#line 310 "dlagts.f"
		if (k >= 3) {
#line 311 "dlagts.f"
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
#line 312 "dlagts.f"
		} else if (k == 2) {
#line 313 "dlagts.f"
		    temp = y[k] - b[k - 1] * y[k - 1];
#line 314 "dlagts.f"
		} else {
#line 315 "dlagts.f"
		    temp = y[k];
#line 316 "dlagts.f"
		}
#line 317 "dlagts.f"
		ak = a[k];
#line 318 "dlagts.f"
		absak = abs(ak);
#line 319 "dlagts.f"
		if (absak < 1.) {
#line 320 "dlagts.f"
		    if (absak < sfmin) {
#line 321 "dlagts.f"
			if (absak == 0. || abs(temp) * sfmin > absak) {
#line 323 "dlagts.f"
			    *info = k;
#line 324 "dlagts.f"
			    return 0;
#line 325 "dlagts.f"
			} else {
#line 326 "dlagts.f"
			    temp *= bignum;
#line 327 "dlagts.f"
			    ak *= bignum;
#line 328 "dlagts.f"
			}
#line 329 "dlagts.f"
		    } else if (abs(temp) > absak * bignum) {
#line 330 "dlagts.f"
			*info = k;
#line 331 "dlagts.f"
			return 0;
#line 332 "dlagts.f"
		    }
#line 333 "dlagts.f"
		}
#line 334 "dlagts.f"
		y[k] = temp / ak;
#line 335 "dlagts.f"
/* L60: */
#line 335 "dlagts.f"
	    }
#line 336 "dlagts.f"
	} else {
#line 337 "dlagts.f"
	    i__1 = *n;
#line 337 "dlagts.f"
	    for (k = 1; k <= i__1; ++k) {
#line 338 "dlagts.f"
		if (k >= 3) {
#line 339 "dlagts.f"
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
#line 340 "dlagts.f"
		} else if (k == 2) {
#line 341 "dlagts.f"
		    temp = y[k] - b[k - 1] * y[k - 1];
#line 342 "dlagts.f"
		} else {
#line 343 "dlagts.f"
		    temp = y[k];
#line 344 "dlagts.f"
		}
#line 345 "dlagts.f"
		ak = a[k];
#line 346 "dlagts.f"
		pert = d_sign(tol, &ak);
#line 347 "dlagts.f"
L70:
#line 348 "dlagts.f"
		absak = abs(ak);
#line 349 "dlagts.f"
		if (absak < 1.) {
#line 350 "dlagts.f"
		    if (absak < sfmin) {
#line 351 "dlagts.f"
			if (absak == 0. || abs(temp) * sfmin > absak) {
#line 353 "dlagts.f"
			    ak += pert;
#line 354 "dlagts.f"
			    pert *= 2;
#line 355 "dlagts.f"
			    goto L70;
#line 356 "dlagts.f"
			} else {
#line 357 "dlagts.f"
			    temp *= bignum;
#line 358 "dlagts.f"
			    ak *= bignum;
#line 359 "dlagts.f"
			}
#line 360 "dlagts.f"
		    } else if (abs(temp) > absak * bignum) {
#line 361 "dlagts.f"
			ak += pert;
#line 362 "dlagts.f"
			pert *= 2;
#line 363 "dlagts.f"
			goto L70;
#line 364 "dlagts.f"
		    }
#line 365 "dlagts.f"
		}
#line 366 "dlagts.f"
		y[k] = temp / ak;
#line 367 "dlagts.f"
/* L80: */
#line 367 "dlagts.f"
	    }
#line 368 "dlagts.f"
	}

#line 370 "dlagts.f"
	for (k = *n; k >= 2; --k) {
#line 371 "dlagts.f"
	    if (in[k - 1] == 0) {
#line 372 "dlagts.f"
		y[k - 1] -= c__[k - 1] * y[k];
#line 373 "dlagts.f"
	    } else {
#line 374 "dlagts.f"
		temp = y[k - 1];
#line 375 "dlagts.f"
		y[k - 1] = y[k];
#line 376 "dlagts.f"
		y[k] = temp - c__[k - 1] * y[k];
#line 377 "dlagts.f"
	    }
#line 378 "dlagts.f"
/* L90: */
#line 378 "dlagts.f"
	}
#line 379 "dlagts.f"
    }

/*     End of DLAGTS */

#line 383 "dlagts.f"
    return 0;
} /* dlagts_ */

