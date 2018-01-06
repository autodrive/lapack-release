#line 1 "slagts.f"
/* slagts.f -- translated by f2c (version 20100827).
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

#line 1 "slagts.f"
/* > \brief \b SLAGTS solves the system of equations (T-λI)x = y or (T-λI)Tx = y,where T is a general tridia
gonal matrix and λ a scalar, using the LU factorization computed by slagtf. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAGTS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagts.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagts.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagts.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, JOB, N */
/*       REAL               TOL */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IN( * ) */
/*       REAL               A( * ), B( * ), C( * ), D( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAGTS may be used to solve one of the systems of equations */
/* > */
/* >    (T - lambda*I)*x = y   or   (T - lambda*I)**T*x = y, */
/* > */
/* > where T is an n by n tridiagonal matrix, for x, following the */
/* > factorization of (T - lambda*I) as */
/* > */
/* >    (T - lambda*I) = P*L*U , */
/* > */
/* > by routine SLAGTF. The choice of equation to be solved is */
/* > controlled by the argument JOB, and in each case there is an option */
/* > to perturb zero or very small diagonal elements of U, this option */
/* > being intended for use in applications such as inverse iteration. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is INTEGER */
/* >          Specifies the job to be performed by SLAGTS as follows: */
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
/* >          A is REAL array, dimension (N) */
/* >          On entry, A must contain the diagonal elements of U as */
/* >          returned from SLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array, dimension (N-1) */
/* >          On entry, B must contain the first super-diagonal elements of */
/* >          U as returned from SLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N-1) */
/* >          On entry, C must contain the sub-diagonal elements of L as */
/* >          returned from SLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N-2) */
/* >          On entry, D must contain the second super-diagonal elements */
/* >          of U as returned from SLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in] IN */
/* > \verbatim */
/* >          IN is INTEGER array, dimension (N) */
/* >          On entry, IN must contain details of the matrix P as returned */
/* >          from SLAGTF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array, dimension (N) */
/* >          On entry, the right hand side vector y. */
/* >          On exit, Y is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[in,out] TOL */
/* > \verbatim */
/* >          TOL is REAL */
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

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slagts_(integer *job, integer *n, doublereal *a, 
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
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

#line 200 "slagts.f"
    /* Parameter adjustments */
#line 200 "slagts.f"
    --y;
#line 200 "slagts.f"
    --in;
#line 200 "slagts.f"
    --d__;
#line 200 "slagts.f"
    --c__;
#line 200 "slagts.f"
    --b;
#line 200 "slagts.f"
    --a;
#line 200 "slagts.f"

#line 200 "slagts.f"
    /* Function Body */
#line 200 "slagts.f"
    *info = 0;
#line 201 "slagts.f"
    if (abs(*job) > 2 || *job == 0) {
#line 202 "slagts.f"
	*info = -1;
#line 203 "slagts.f"
    } else if (*n < 0) {
#line 204 "slagts.f"
	*info = -2;
#line 205 "slagts.f"
    }
#line 206 "slagts.f"
    if (*info != 0) {
#line 207 "slagts.f"
	i__1 = -(*info);
#line 207 "slagts.f"
	xerbla_("SLAGTS", &i__1, (ftnlen)6);
#line 208 "slagts.f"
	return 0;
#line 209 "slagts.f"
    }

#line 211 "slagts.f"
    if (*n == 0) {
#line 211 "slagts.f"
	return 0;
#line 211 "slagts.f"
    }

#line 214 "slagts.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 215 "slagts.f"
    sfmin = slamch_("Safe minimum", (ftnlen)12);
#line 216 "slagts.f"
    bignum = 1. / sfmin;

#line 218 "slagts.f"
    if (*job < 0) {
#line 219 "slagts.f"
	if (*tol <= 0.) {
#line 220 "slagts.f"
	    *tol = abs(a[1]);
#line 221 "slagts.f"
	    if (*n > 1) {
/* Computing MAX */
#line 221 "slagts.f"
		d__1 = *tol, d__2 = abs(a[2]), d__1 = max(d__1,d__2), d__2 = 
			abs(b[1]);
#line 221 "slagts.f"
		*tol = max(d__1,d__2);
#line 221 "slagts.f"
	    }
#line 223 "slagts.f"
	    i__1 = *n;
#line 223 "slagts.f"
	    for (k = 3; k <= i__1; ++k) {
/* Computing MAX */
#line 224 "slagts.f"
		d__4 = *tol, d__5 = (d__1 = a[k], abs(d__1)), d__4 = max(d__4,
			d__5), d__5 = (d__2 = b[k - 1], abs(d__2)), d__4 = 
			max(d__4,d__5), d__5 = (d__3 = d__[k - 2], abs(d__3));
#line 224 "slagts.f"
		*tol = max(d__4,d__5);
#line 226 "slagts.f"
/* L10: */
#line 226 "slagts.f"
	    }
#line 227 "slagts.f"
	    *tol *= eps;
#line 228 "slagts.f"
	    if (*tol == 0.) {
#line 228 "slagts.f"
		*tol = eps;
#line 228 "slagts.f"
	    }
#line 230 "slagts.f"
	}
#line 231 "slagts.f"
    }

#line 233 "slagts.f"
    if (abs(*job) == 1) {
#line 234 "slagts.f"
	i__1 = *n;
#line 234 "slagts.f"
	for (k = 2; k <= i__1; ++k) {
#line 235 "slagts.f"
	    if (in[k - 1] == 0) {
#line 236 "slagts.f"
		y[k] -= c__[k - 1] * y[k - 1];
#line 237 "slagts.f"
	    } else {
#line 238 "slagts.f"
		temp = y[k - 1];
#line 239 "slagts.f"
		y[k - 1] = y[k];
#line 240 "slagts.f"
		y[k] = temp - c__[k - 1] * y[k];
#line 241 "slagts.f"
	    }
#line 242 "slagts.f"
/* L20: */
#line 242 "slagts.f"
	}
#line 243 "slagts.f"
	if (*job == 1) {
#line 244 "slagts.f"
	    for (k = *n; k >= 1; --k) {
#line 245 "slagts.f"
		if (k <= *n - 2) {
#line 246 "slagts.f"
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
#line 247 "slagts.f"
		} else if (k == *n - 1) {
#line 248 "slagts.f"
		    temp = y[k] - b[k] * y[k + 1];
#line 249 "slagts.f"
		} else {
#line 250 "slagts.f"
		    temp = y[k];
#line 251 "slagts.f"
		}
#line 252 "slagts.f"
		ak = a[k];
#line 253 "slagts.f"
		absak = abs(ak);
#line 254 "slagts.f"
		if (absak < 1.) {
#line 255 "slagts.f"
		    if (absak < sfmin) {
#line 256 "slagts.f"
			if (absak == 0. || abs(temp) * sfmin > absak) {
#line 258 "slagts.f"
			    *info = k;
#line 259 "slagts.f"
			    return 0;
#line 260 "slagts.f"
			} else {
#line 261 "slagts.f"
			    temp *= bignum;
#line 262 "slagts.f"
			    ak *= bignum;
#line 263 "slagts.f"
			}
#line 264 "slagts.f"
		    } else if (abs(temp) > absak * bignum) {
#line 265 "slagts.f"
			*info = k;
#line 266 "slagts.f"
			return 0;
#line 267 "slagts.f"
		    }
#line 268 "slagts.f"
		}
#line 269 "slagts.f"
		y[k] = temp / ak;
#line 270 "slagts.f"
/* L30: */
#line 270 "slagts.f"
	    }
#line 271 "slagts.f"
	} else {
#line 272 "slagts.f"
	    for (k = *n; k >= 1; --k) {
#line 273 "slagts.f"
		if (k <= *n - 2) {
#line 274 "slagts.f"
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
#line 275 "slagts.f"
		} else if (k == *n - 1) {
#line 276 "slagts.f"
		    temp = y[k] - b[k] * y[k + 1];
#line 277 "slagts.f"
		} else {
#line 278 "slagts.f"
		    temp = y[k];
#line 279 "slagts.f"
		}
#line 280 "slagts.f"
		ak = a[k];
#line 281 "slagts.f"
		pert = d_sign(tol, &ak);
#line 282 "slagts.f"
L40:
#line 283 "slagts.f"
		absak = abs(ak);
#line 284 "slagts.f"
		if (absak < 1.) {
#line 285 "slagts.f"
		    if (absak < sfmin) {
#line 286 "slagts.f"
			if (absak == 0. || abs(temp) * sfmin > absak) {
#line 288 "slagts.f"
			    ak += pert;
#line 289 "slagts.f"
			    pert *= 2;
#line 290 "slagts.f"
			    goto L40;
#line 291 "slagts.f"
			} else {
#line 292 "slagts.f"
			    temp *= bignum;
#line 293 "slagts.f"
			    ak *= bignum;
#line 294 "slagts.f"
			}
#line 295 "slagts.f"
		    } else if (abs(temp) > absak * bignum) {
#line 296 "slagts.f"
			ak += pert;
#line 297 "slagts.f"
			pert *= 2;
#line 298 "slagts.f"
			goto L40;
#line 299 "slagts.f"
		    }
#line 300 "slagts.f"
		}
#line 301 "slagts.f"
		y[k] = temp / ak;
#line 302 "slagts.f"
/* L50: */
#line 302 "slagts.f"
	    }
#line 303 "slagts.f"
	}
#line 304 "slagts.f"
    } else {

/*        Come to here if  JOB = 2 or -2 */

#line 308 "slagts.f"
	if (*job == 2) {
#line 309 "slagts.f"
	    i__1 = *n;
#line 309 "slagts.f"
	    for (k = 1; k <= i__1; ++k) {
#line 310 "slagts.f"
		if (k >= 3) {
#line 311 "slagts.f"
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
#line 312 "slagts.f"
		} else if (k == 2) {
#line 313 "slagts.f"
		    temp = y[k] - b[k - 1] * y[k - 1];
#line 314 "slagts.f"
		} else {
#line 315 "slagts.f"
		    temp = y[k];
#line 316 "slagts.f"
		}
#line 317 "slagts.f"
		ak = a[k];
#line 318 "slagts.f"
		absak = abs(ak);
#line 319 "slagts.f"
		if (absak < 1.) {
#line 320 "slagts.f"
		    if (absak < sfmin) {
#line 321 "slagts.f"
			if (absak == 0. || abs(temp) * sfmin > absak) {
#line 323 "slagts.f"
			    *info = k;
#line 324 "slagts.f"
			    return 0;
#line 325 "slagts.f"
			} else {
#line 326 "slagts.f"
			    temp *= bignum;
#line 327 "slagts.f"
			    ak *= bignum;
#line 328 "slagts.f"
			}
#line 329 "slagts.f"
		    } else if (abs(temp) > absak * bignum) {
#line 330 "slagts.f"
			*info = k;
#line 331 "slagts.f"
			return 0;
#line 332 "slagts.f"
		    }
#line 333 "slagts.f"
		}
#line 334 "slagts.f"
		y[k] = temp / ak;
#line 335 "slagts.f"
/* L60: */
#line 335 "slagts.f"
	    }
#line 336 "slagts.f"
	} else {
#line 337 "slagts.f"
	    i__1 = *n;
#line 337 "slagts.f"
	    for (k = 1; k <= i__1; ++k) {
#line 338 "slagts.f"
		if (k >= 3) {
#line 339 "slagts.f"
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
#line 340 "slagts.f"
		} else if (k == 2) {
#line 341 "slagts.f"
		    temp = y[k] - b[k - 1] * y[k - 1];
#line 342 "slagts.f"
		} else {
#line 343 "slagts.f"
		    temp = y[k];
#line 344 "slagts.f"
		}
#line 345 "slagts.f"
		ak = a[k];
#line 346 "slagts.f"
		pert = d_sign(tol, &ak);
#line 347 "slagts.f"
L70:
#line 348 "slagts.f"
		absak = abs(ak);
#line 349 "slagts.f"
		if (absak < 1.) {
#line 350 "slagts.f"
		    if (absak < sfmin) {
#line 351 "slagts.f"
			if (absak == 0. || abs(temp) * sfmin > absak) {
#line 353 "slagts.f"
			    ak += pert;
#line 354 "slagts.f"
			    pert *= 2;
#line 355 "slagts.f"
			    goto L70;
#line 356 "slagts.f"
			} else {
#line 357 "slagts.f"
			    temp *= bignum;
#line 358 "slagts.f"
			    ak *= bignum;
#line 359 "slagts.f"
			}
#line 360 "slagts.f"
		    } else if (abs(temp) > absak * bignum) {
#line 361 "slagts.f"
			ak += pert;
#line 362 "slagts.f"
			pert *= 2;
#line 363 "slagts.f"
			goto L70;
#line 364 "slagts.f"
		    }
#line 365 "slagts.f"
		}
#line 366 "slagts.f"
		y[k] = temp / ak;
#line 367 "slagts.f"
/* L80: */
#line 367 "slagts.f"
	    }
#line 368 "slagts.f"
	}

#line 370 "slagts.f"
	for (k = *n; k >= 2; --k) {
#line 371 "slagts.f"
	    if (in[k - 1] == 0) {
#line 372 "slagts.f"
		y[k - 1] -= c__[k - 1] * y[k];
#line 373 "slagts.f"
	    } else {
#line 374 "slagts.f"
		temp = y[k - 1];
#line 375 "slagts.f"
		y[k - 1] = y[k];
#line 376 "slagts.f"
		y[k] = temp - c__[k - 1] * y[k];
#line 377 "slagts.f"
	    }
#line 378 "slagts.f"
/* L90: */
#line 378 "slagts.f"
	}
#line 379 "slagts.f"
    }

/*     End of SLAGTS */

#line 383 "slagts.f"
    return 0;
} /* slagts_ */

