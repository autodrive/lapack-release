#line 1 "slaed6.f"
/* slaed6.f -- translated by f2c (version 20100827).
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

#line 1 "slaed6.f"
/* > \brief \b SLAED6 used by sstedc. Computes one Newton step in solution of the secular equation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAED6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed6.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed6.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed6.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            ORGATI */
/*       INTEGER            INFO, KNITER */
/*       REAL               FINIT, RHO, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( 3 ), Z( 3 ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAED6 computes the positive or negative root (closest to the origin) */
/* > of */
/* >                  z(1)        z(2)        z(3) */
/* > f(x) =   rho + --------- + ---------- + --------- */
/* >                 d(1)-x      d(2)-x      d(3)-x */
/* > */
/* > It is assumed that */
/* > */
/* >       if ORGATI = .true. the root is between d(2) and d(3); */
/* >       otherwise it is between d(1) and d(2) */
/* > */
/* > This routine will be called by SLAED4 when necessary. In most cases, */
/* > the root sought is the smallest in magnitude, though it might not be */
/* > in some extremely rare situations. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] KNITER */
/* > \verbatim */
/* >          KNITER is INTEGER */
/* >               Refer to SLAED4 for its significance. */
/* > \endverbatim */
/* > */
/* > \param[in] ORGATI */
/* > \verbatim */
/* >          ORGATI is LOGICAL */
/* >               If ORGATI is true, the needed root is between d(2) and */
/* >               d(3); otherwise it is between d(1) and d(2).  See */
/* >               SLAED4 for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is REAL */
/* >               Refer to the equation f(x) above. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (3) */
/* >               D satisfies d(1) < d(2) < d(3). */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (3) */
/* >               Each of the elements in z must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in] FINIT */
/* > \verbatim */
/* >          FINIT is REAL */
/* >               The value of f at 0. It is more accurate than the one */
/* >               evaluated inside this routine (if someone wants to do */
/* >               so). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL */
/* >               The root of the equation f(x). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >               = 0: successful exit */
/* >               > 0: if INFO = 1, failure to converge */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  10/02/03: This version has a few statements commented out for thread */
/* >  safety (machine parameters are computed on each entry). SJH. */
/* > */
/* >  05/10/06: Modified from a new version of Ren-Cang Li, use */
/* >     Gragg-Thornton-Warner cubic convergent scheme for better stability. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ren-Cang Li, Computer Science Division, University of California */
/* >     at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slaed6_(integer *kniter, logical *orgati, doublereal *
	rho, doublereal *d__, doublereal *z__, doublereal *finit, doublereal *
	tau, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal a, b, c__, f;
    static integer i__;
    static doublereal fc, df, ddf, lbd, eta, ubd, eps, base;
    static integer iter;
    static doublereal temp, temp1, temp2, temp3, temp4;
    static logical scale;
    static integer niter;
    static doublereal small1, small2, sminv1, sminv2, dscale[3], sclfac;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal zscale[3], erretm, sclinv;


/*  -- LAPACK computational routine (version 3.4.2) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 186 "slaed6.f"
    /* Parameter adjustments */
#line 186 "slaed6.f"
    --z__;
#line 186 "slaed6.f"
    --d__;
#line 186 "slaed6.f"

#line 186 "slaed6.f"
    /* Function Body */
#line 186 "slaed6.f"
    *info = 0;

#line 188 "slaed6.f"
    if (*orgati) {
#line 189 "slaed6.f"
	lbd = d__[2];
#line 190 "slaed6.f"
	ubd = d__[3];
#line 191 "slaed6.f"
    } else {
#line 192 "slaed6.f"
	lbd = d__[1];
#line 193 "slaed6.f"
	ubd = d__[2];
#line 194 "slaed6.f"
    }
#line 195 "slaed6.f"
    if (*finit < 0.) {
#line 196 "slaed6.f"
	lbd = 0.;
#line 197 "slaed6.f"
    } else {
#line 198 "slaed6.f"
	ubd = 0.;
#line 199 "slaed6.f"
    }

#line 201 "slaed6.f"
    niter = 1;
#line 202 "slaed6.f"
    *tau = 0.;
#line 203 "slaed6.f"
    if (*kniter == 2) {
#line 204 "slaed6.f"
	if (*orgati) {
#line 205 "slaed6.f"
	    temp = (d__[3] - d__[2]) / 2.;
#line 206 "slaed6.f"
	    c__ = *rho + z__[1] / (d__[1] - d__[2] - temp);
#line 207 "slaed6.f"
	    a = c__ * (d__[2] + d__[3]) + z__[2] + z__[3];
#line 208 "slaed6.f"
	    b = c__ * d__[2] * d__[3] + z__[2] * d__[3] + z__[3] * d__[2];
#line 209 "slaed6.f"
	} else {
#line 210 "slaed6.f"
	    temp = (d__[1] - d__[2]) / 2.;
#line 211 "slaed6.f"
	    c__ = *rho + z__[3] / (d__[3] - d__[2] - temp);
#line 212 "slaed6.f"
	    a = c__ * (d__[1] + d__[2]) + z__[1] + z__[2];
#line 213 "slaed6.f"
	    b = c__ * d__[1] * d__[2] + z__[1] * d__[2] + z__[2] * d__[1];
#line 214 "slaed6.f"
	}
/* Computing MAX */
#line 215 "slaed6.f"
	d__1 = abs(a), d__2 = abs(b), d__1 = max(d__1,d__2), d__2 = abs(c__);
#line 215 "slaed6.f"
	temp = max(d__1,d__2);
#line 216 "slaed6.f"
	a /= temp;
#line 217 "slaed6.f"
	b /= temp;
#line 218 "slaed6.f"
	c__ /= temp;
#line 219 "slaed6.f"
	if (c__ == 0.) {
#line 220 "slaed6.f"
	    *tau = b / a;
#line 221 "slaed6.f"
	} else if (a <= 0.) {
#line 222 "slaed6.f"
	    *tau = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
		    c__ * 2.);
#line 223 "slaed6.f"
	} else {
#line 224 "slaed6.f"
	    *tau = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))
		    ));
#line 225 "slaed6.f"
	}
#line 226 "slaed6.f"
	if (*tau < lbd || *tau > ubd) {
#line 226 "slaed6.f"
	    *tau = (lbd + ubd) / 2.;
#line 226 "slaed6.f"
	}
#line 228 "slaed6.f"
	if (d__[1] == *tau || d__[2] == *tau || d__[3] == *tau) {
#line 229 "slaed6.f"
	    *tau = 0.;
#line 230 "slaed6.f"
	} else {
#line 231 "slaed6.f"
	    temp = *finit + *tau * z__[1] / (d__[1] * (d__[1] - *tau)) + *tau 
		    * z__[2] / (d__[2] * (d__[2] - *tau)) + *tau * z__[3] / (
		    d__[3] * (d__[3] - *tau));
#line 234 "slaed6.f"
	    if (temp <= 0.) {
#line 235 "slaed6.f"
		lbd = *tau;
#line 236 "slaed6.f"
	    } else {
#line 237 "slaed6.f"
		ubd = *tau;
#line 238 "slaed6.f"
	    }
#line 239 "slaed6.f"
	    if (abs(*finit) <= abs(temp)) {
#line 239 "slaed6.f"
		*tau = 0.;
#line 239 "slaed6.f"
	    }
#line 241 "slaed6.f"
	}
#line 242 "slaed6.f"
    }

/*     get machine parameters for possible scaling to avoid overflow */

/*     modified by Sven: parameters SMALL1, SMINV1, SMALL2, */
/*     SMINV2, EPS are not SAVEd anymore between one call to the */
/*     others but recomputed at each call */

#line 250 "slaed6.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 251 "slaed6.f"
    base = slamch_("Base", (ftnlen)4);
#line 252 "slaed6.f"
    i__1 = (integer) (log(slamch_("SafMin", (ftnlen)6)) / log(base) / 3.);
#line 252 "slaed6.f"
    small1 = pow_di(&base, &i__1);
#line 254 "slaed6.f"
    sminv1 = 1. / small1;
#line 255 "slaed6.f"
    small2 = small1 * small1;
#line 256 "slaed6.f"
    sminv2 = sminv1 * sminv1;

/*     Determine if scaling of inputs necessary to avoid overflow */
/*     when computing 1/TEMP**3 */

#line 261 "slaed6.f"
    if (*orgati) {
/* Computing MIN */
#line 262 "slaed6.f"
	d__3 = (d__1 = d__[2] - *tau, abs(d__1)), d__4 = (d__2 = d__[3] - *
		tau, abs(d__2));
#line 262 "slaed6.f"
	temp = min(d__3,d__4);
#line 263 "slaed6.f"
    } else {
/* Computing MIN */
#line 264 "slaed6.f"
	d__3 = (d__1 = d__[1] - *tau, abs(d__1)), d__4 = (d__2 = d__[2] - *
		tau, abs(d__2));
#line 264 "slaed6.f"
	temp = min(d__3,d__4);
#line 265 "slaed6.f"
    }
#line 266 "slaed6.f"
    scale = FALSE_;
#line 267 "slaed6.f"
    if (temp <= small1) {
#line 268 "slaed6.f"
	scale = TRUE_;
#line 269 "slaed6.f"
	if (temp <= small2) {

/*        Scale up by power of radix nearest 1/SAFMIN**(2/3) */

#line 273 "slaed6.f"
	    sclfac = sminv2;
#line 274 "slaed6.f"
	    sclinv = small2;
#line 275 "slaed6.f"
	} else {

/*        Scale up by power of radix nearest 1/SAFMIN**(1/3) */

#line 279 "slaed6.f"
	    sclfac = sminv1;
#line 280 "slaed6.f"
	    sclinv = small1;
#line 281 "slaed6.f"
	}

/*        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1) */

#line 285 "slaed6.f"
	for (i__ = 1; i__ <= 3; ++i__) {
#line 286 "slaed6.f"
	    dscale[i__ - 1] = d__[i__] * sclfac;
#line 287 "slaed6.f"
	    zscale[i__ - 1] = z__[i__] * sclfac;
#line 288 "slaed6.f"
/* L10: */
#line 288 "slaed6.f"
	}
#line 289 "slaed6.f"
	*tau *= sclfac;
#line 290 "slaed6.f"
	lbd *= sclfac;
#line 291 "slaed6.f"
	ubd *= sclfac;
#line 292 "slaed6.f"
    } else {

/*        Copy D and Z to DSCALE and ZSCALE */

#line 296 "slaed6.f"
	for (i__ = 1; i__ <= 3; ++i__) {
#line 297 "slaed6.f"
	    dscale[i__ - 1] = d__[i__];
#line 298 "slaed6.f"
	    zscale[i__ - 1] = z__[i__];
#line 299 "slaed6.f"
/* L20: */
#line 299 "slaed6.f"
	}
#line 300 "slaed6.f"
    }

#line 302 "slaed6.f"
    fc = 0.;
#line 303 "slaed6.f"
    df = 0.;
#line 304 "slaed6.f"
    ddf = 0.;
#line 305 "slaed6.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 306 "slaed6.f"
	temp = 1. / (dscale[i__ - 1] - *tau);
#line 307 "slaed6.f"
	temp1 = zscale[i__ - 1] * temp;
#line 308 "slaed6.f"
	temp2 = temp1 * temp;
#line 309 "slaed6.f"
	temp3 = temp2 * temp;
#line 310 "slaed6.f"
	fc += temp1 / dscale[i__ - 1];
#line 311 "slaed6.f"
	df += temp2;
#line 312 "slaed6.f"
	ddf += temp3;
#line 313 "slaed6.f"
/* L30: */
#line 313 "slaed6.f"
    }
#line 314 "slaed6.f"
    f = *finit + *tau * fc;

#line 316 "slaed6.f"
    if (abs(f) <= 0.) {
#line 316 "slaed6.f"
	goto L60;
#line 316 "slaed6.f"
    }
#line 318 "slaed6.f"
    if (f <= 0.) {
#line 319 "slaed6.f"
	lbd = *tau;
#line 320 "slaed6.f"
    } else {
#line 321 "slaed6.f"
	ubd = *tau;
#line 322 "slaed6.f"
    }

/*        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent */
/*                            scheme */

/*     It is not hard to see that */

/*           1) Iterations will go up monotonically */
/*              if FINIT < 0; */

/*           2) Iterations will go down monotonically */
/*              if FINIT > 0. */

#line 335 "slaed6.f"
    iter = niter + 1;

#line 337 "slaed6.f"
    for (niter = iter; niter <= 40; ++niter) {

#line 339 "slaed6.f"
	if (*orgati) {
#line 340 "slaed6.f"
	    temp1 = dscale[1] - *tau;
#line 341 "slaed6.f"
	    temp2 = dscale[2] - *tau;
#line 342 "slaed6.f"
	} else {
#line 343 "slaed6.f"
	    temp1 = dscale[0] - *tau;
#line 344 "slaed6.f"
	    temp2 = dscale[1] - *tau;
#line 345 "slaed6.f"
	}
#line 346 "slaed6.f"
	a = (temp1 + temp2) * f - temp1 * temp2 * df;
#line 347 "slaed6.f"
	b = temp1 * temp2 * f;
#line 348 "slaed6.f"
	c__ = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
/* Computing MAX */
#line 349 "slaed6.f"
	d__1 = abs(a), d__2 = abs(b), d__1 = max(d__1,d__2), d__2 = abs(c__);
#line 349 "slaed6.f"
	temp = max(d__1,d__2);
#line 350 "slaed6.f"
	a /= temp;
#line 351 "slaed6.f"
	b /= temp;
#line 352 "slaed6.f"
	c__ /= temp;
#line 353 "slaed6.f"
	if (c__ == 0.) {
#line 354 "slaed6.f"
	    eta = b / a;
#line 355 "slaed6.f"
	} else if (a <= 0.) {
#line 356 "slaed6.f"
	    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ 
		    * 2.);
#line 357 "slaed6.f"
	} else {
#line 358 "slaed6.f"
	    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))
		    );
#line 359 "slaed6.f"
	}
#line 360 "slaed6.f"
	if (f * eta >= 0.) {
#line 361 "slaed6.f"
	    eta = -f / df;
#line 362 "slaed6.f"
	}

#line 364 "slaed6.f"
	*tau += eta;
#line 365 "slaed6.f"
	if (*tau < lbd || *tau > ubd) {
#line 365 "slaed6.f"
	    *tau = (lbd + ubd) / 2.;
#line 365 "slaed6.f"
	}

#line 368 "slaed6.f"
	fc = 0.;
#line 369 "slaed6.f"
	erretm = 0.;
#line 370 "slaed6.f"
	df = 0.;
#line 371 "slaed6.f"
	ddf = 0.;
#line 372 "slaed6.f"
	for (i__ = 1; i__ <= 3; ++i__) {
#line 373 "slaed6.f"
	    if (dscale[i__ - 1] - *tau != 0.) {
#line 374 "slaed6.f"
		temp = 1. / (dscale[i__ - 1] - *tau);
#line 375 "slaed6.f"
		temp1 = zscale[i__ - 1] * temp;
#line 376 "slaed6.f"
		temp2 = temp1 * temp;
#line 377 "slaed6.f"
		temp3 = temp2 * temp;
#line 378 "slaed6.f"
		temp4 = temp1 / dscale[i__ - 1];
#line 379 "slaed6.f"
		fc += temp4;
#line 380 "slaed6.f"
		erretm += abs(temp4);
#line 381 "slaed6.f"
		df += temp2;
#line 382 "slaed6.f"
		ddf += temp3;
#line 383 "slaed6.f"
	    } else {
#line 384 "slaed6.f"
		goto L60;
#line 385 "slaed6.f"
	    }
#line 386 "slaed6.f"
/* L40: */
#line 386 "slaed6.f"
	}
#line 387 "slaed6.f"
	f = *finit + *tau * fc;
#line 388 "slaed6.f"
	erretm = (abs(*finit) + abs(*tau) * erretm) * 8. + abs(*tau) * df;
#line 390 "slaed6.f"
	if (abs(f) <= eps * erretm) {
#line 390 "slaed6.f"
	    goto L60;
#line 390 "slaed6.f"
	}
#line 392 "slaed6.f"
	if (f <= 0.) {
#line 393 "slaed6.f"
	    lbd = *tau;
#line 394 "slaed6.f"
	} else {
#line 395 "slaed6.f"
	    ubd = *tau;
#line 396 "slaed6.f"
	}
#line 397 "slaed6.f"
/* L50: */
#line 397 "slaed6.f"
    }
#line 398 "slaed6.f"
    *info = 1;
#line 399 "slaed6.f"
L60:

/*     Undo scaling */

#line 403 "slaed6.f"
    if (scale) {
#line 403 "slaed6.f"
	*tau *= sclinv;
#line 403 "slaed6.f"
    }
#line 405 "slaed6.f"
    return 0;

/*     End of SLAED6 */

} /* slaed6_ */

