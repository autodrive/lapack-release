#line 1 "dlaed6.f"
/* dlaed6.f -- translated by f2c (version 20100827).
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

#line 1 "dlaed6.f"
/* > \brief \b DLAED6 used by sstedc. Computes one Newton step in solution of the secular equation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAED6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed6.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed6.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed6.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            ORGATI */
/*       INTEGER            INFO, KNITER */
/*       DOUBLE PRECISION   FINIT, RHO, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( 3 ), Z( 3 ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAED6 computes the positive or negative root (closest to the origin) */
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
/* > This routine will be called by DLAED4 when necessary. In most cases, */
/* > the root sought is the smallest in magnitude, though it might not be */
/* > in some extremely rare situations. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] KNITER */
/* > \verbatim */
/* >          KNITER is INTEGER */
/* >               Refer to DLAED4 for its significance. */
/* > \endverbatim */
/* > */
/* > \param[in] ORGATI */
/* > \verbatim */
/* >          ORGATI is LOGICAL */
/* >               If ORGATI is true, the needed root is between d(2) and */
/* >               d(3); otherwise it is between d(1) and d(2).  See */
/* >               DLAED4 for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >               Refer to the equation f(x) above. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (3) */
/* >               D satisfies d(1) < d(2) < d(3). */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (3) */
/* >               Each of the elements in z must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in] FINIT */
/* > \verbatim */
/* >          FINIT is DOUBLE PRECISION */
/* >               The value of f at 0. It is more accurate than the one */
/* >               evaluated inside this routine (if someone wants to do */
/* >               so). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION */
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

/* > \date December 2016 */

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
/* Subroutine */ int dlaed6_(integer *kniter, logical *orgati, doublereal *
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
    static doublereal small1, small2, sminv1, sminv2;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal dscale[3], sclfac, zscale[3], erretm, sclinv;


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
/*     .. External Functions .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 186 "dlaed6.f"
    /* Parameter adjustments */
#line 186 "dlaed6.f"
    --z__;
#line 186 "dlaed6.f"
    --d__;
#line 186 "dlaed6.f"

#line 186 "dlaed6.f"
    /* Function Body */
#line 186 "dlaed6.f"
    *info = 0;

#line 188 "dlaed6.f"
    if (*orgati) {
#line 189 "dlaed6.f"
	lbd = d__[2];
#line 190 "dlaed6.f"
	ubd = d__[3];
#line 191 "dlaed6.f"
    } else {
#line 192 "dlaed6.f"
	lbd = d__[1];
#line 193 "dlaed6.f"
	ubd = d__[2];
#line 194 "dlaed6.f"
    }
#line 195 "dlaed6.f"
    if (*finit < 0.) {
#line 196 "dlaed6.f"
	lbd = 0.;
#line 197 "dlaed6.f"
    } else {
#line 198 "dlaed6.f"
	ubd = 0.;
#line 199 "dlaed6.f"
    }

#line 201 "dlaed6.f"
    niter = 1;
#line 202 "dlaed6.f"
    *tau = 0.;
#line 203 "dlaed6.f"
    if (*kniter == 2) {
#line 204 "dlaed6.f"
	if (*orgati) {
#line 205 "dlaed6.f"
	    temp = (d__[3] - d__[2]) / 2.;
#line 206 "dlaed6.f"
	    c__ = *rho + z__[1] / (d__[1] - d__[2] - temp);
#line 207 "dlaed6.f"
	    a = c__ * (d__[2] + d__[3]) + z__[2] + z__[3];
#line 208 "dlaed6.f"
	    b = c__ * d__[2] * d__[3] + z__[2] * d__[3] + z__[3] * d__[2];
#line 209 "dlaed6.f"
	} else {
#line 210 "dlaed6.f"
	    temp = (d__[1] - d__[2]) / 2.;
#line 211 "dlaed6.f"
	    c__ = *rho + z__[3] / (d__[3] - d__[2] - temp);
#line 212 "dlaed6.f"
	    a = c__ * (d__[1] + d__[2]) + z__[1] + z__[2];
#line 213 "dlaed6.f"
	    b = c__ * d__[1] * d__[2] + z__[1] * d__[2] + z__[2] * d__[1];
#line 214 "dlaed6.f"
	}
/* Computing MAX */
#line 215 "dlaed6.f"
	d__1 = abs(a), d__2 = abs(b), d__1 = max(d__1,d__2), d__2 = abs(c__);
#line 215 "dlaed6.f"
	temp = max(d__1,d__2);
#line 216 "dlaed6.f"
	a /= temp;
#line 217 "dlaed6.f"
	b /= temp;
#line 218 "dlaed6.f"
	c__ /= temp;
#line 219 "dlaed6.f"
	if (c__ == 0.) {
#line 220 "dlaed6.f"
	    *tau = b / a;
#line 221 "dlaed6.f"
	} else if (a <= 0.) {
#line 222 "dlaed6.f"
	    *tau = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
		    c__ * 2.);
#line 223 "dlaed6.f"
	} else {
#line 224 "dlaed6.f"
	    *tau = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))
		    ));
#line 225 "dlaed6.f"
	}
#line 226 "dlaed6.f"
	if (*tau < lbd || *tau > ubd) {
#line 226 "dlaed6.f"
	    *tau = (lbd + ubd) / 2.;
#line 226 "dlaed6.f"
	}
#line 228 "dlaed6.f"
	if (d__[1] == *tau || d__[2] == *tau || d__[3] == *tau) {
#line 229 "dlaed6.f"
	    *tau = 0.;
#line 230 "dlaed6.f"
	} else {
#line 231 "dlaed6.f"
	    temp = *finit + *tau * z__[1] / (d__[1] * (d__[1] - *tau)) + *tau 
		    * z__[2] / (d__[2] * (d__[2] - *tau)) + *tau * z__[3] / (
		    d__[3] * (d__[3] - *tau));
#line 234 "dlaed6.f"
	    if (temp <= 0.) {
#line 235 "dlaed6.f"
		lbd = *tau;
#line 236 "dlaed6.f"
	    } else {
#line 237 "dlaed6.f"
		ubd = *tau;
#line 238 "dlaed6.f"
	    }
#line 239 "dlaed6.f"
	    if (abs(*finit) <= abs(temp)) {
#line 239 "dlaed6.f"
		*tau = 0.;
#line 239 "dlaed6.f"
	    }
#line 241 "dlaed6.f"
	}
#line 242 "dlaed6.f"
    }

/*     get machine parameters for possible scaling to avoid overflow */

/*     modified by Sven: parameters SMALL1, SMINV1, SMALL2, */
/*     SMINV2, EPS are not SAVEd anymore between one call to the */
/*     others but recomputed at each call */

#line 250 "dlaed6.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 251 "dlaed6.f"
    base = dlamch_("Base", (ftnlen)4);
#line 252 "dlaed6.f"
    i__1 = (integer) (log(dlamch_("SafMin", (ftnlen)6)) / log(base) / 3.);
#line 252 "dlaed6.f"
    small1 = pow_di(&base, &i__1);
#line 254 "dlaed6.f"
    sminv1 = 1. / small1;
#line 255 "dlaed6.f"
    small2 = small1 * small1;
#line 256 "dlaed6.f"
    sminv2 = sminv1 * sminv1;

/*     Determine if scaling of inputs necessary to avoid overflow */
/*     when computing 1/TEMP**3 */

#line 261 "dlaed6.f"
    if (*orgati) {
/* Computing MIN */
#line 262 "dlaed6.f"
	d__3 = (d__1 = d__[2] - *tau, abs(d__1)), d__4 = (d__2 = d__[3] - *
		tau, abs(d__2));
#line 262 "dlaed6.f"
	temp = min(d__3,d__4);
#line 263 "dlaed6.f"
    } else {
/* Computing MIN */
#line 264 "dlaed6.f"
	d__3 = (d__1 = d__[1] - *tau, abs(d__1)), d__4 = (d__2 = d__[2] - *
		tau, abs(d__2));
#line 264 "dlaed6.f"
	temp = min(d__3,d__4);
#line 265 "dlaed6.f"
    }
#line 266 "dlaed6.f"
    scale = FALSE_;
#line 267 "dlaed6.f"
    if (temp <= small1) {
#line 268 "dlaed6.f"
	scale = TRUE_;
#line 269 "dlaed6.f"
	if (temp <= small2) {

/*        Scale up by power of radix nearest 1/SAFMIN**(2/3) */

#line 273 "dlaed6.f"
	    sclfac = sminv2;
#line 274 "dlaed6.f"
	    sclinv = small2;
#line 275 "dlaed6.f"
	} else {

/*        Scale up by power of radix nearest 1/SAFMIN**(1/3) */

#line 279 "dlaed6.f"
	    sclfac = sminv1;
#line 280 "dlaed6.f"
	    sclinv = small1;
#line 281 "dlaed6.f"
	}

/*        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1) */

#line 285 "dlaed6.f"
	for (i__ = 1; i__ <= 3; ++i__) {
#line 286 "dlaed6.f"
	    dscale[i__ - 1] = d__[i__] * sclfac;
#line 287 "dlaed6.f"
	    zscale[i__ - 1] = z__[i__] * sclfac;
#line 288 "dlaed6.f"
/* L10: */
#line 288 "dlaed6.f"
	}
#line 289 "dlaed6.f"
	*tau *= sclfac;
#line 290 "dlaed6.f"
	lbd *= sclfac;
#line 291 "dlaed6.f"
	ubd *= sclfac;
#line 292 "dlaed6.f"
    } else {

/*        Copy D and Z to DSCALE and ZSCALE */

#line 296 "dlaed6.f"
	for (i__ = 1; i__ <= 3; ++i__) {
#line 297 "dlaed6.f"
	    dscale[i__ - 1] = d__[i__];
#line 298 "dlaed6.f"
	    zscale[i__ - 1] = z__[i__];
#line 299 "dlaed6.f"
/* L20: */
#line 299 "dlaed6.f"
	}
#line 300 "dlaed6.f"
    }

#line 302 "dlaed6.f"
    fc = 0.;
#line 303 "dlaed6.f"
    df = 0.;
#line 304 "dlaed6.f"
    ddf = 0.;
#line 305 "dlaed6.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 306 "dlaed6.f"
	temp = 1. / (dscale[i__ - 1] - *tau);
#line 307 "dlaed6.f"
	temp1 = zscale[i__ - 1] * temp;
#line 308 "dlaed6.f"
	temp2 = temp1 * temp;
#line 309 "dlaed6.f"
	temp3 = temp2 * temp;
#line 310 "dlaed6.f"
	fc += temp1 / dscale[i__ - 1];
#line 311 "dlaed6.f"
	df += temp2;
#line 312 "dlaed6.f"
	ddf += temp3;
#line 313 "dlaed6.f"
/* L30: */
#line 313 "dlaed6.f"
    }
#line 314 "dlaed6.f"
    f = *finit + *tau * fc;

#line 316 "dlaed6.f"
    if (abs(f) <= 0.) {
#line 316 "dlaed6.f"
	goto L60;
#line 316 "dlaed6.f"
    }
#line 318 "dlaed6.f"
    if (f <= 0.) {
#line 319 "dlaed6.f"
	lbd = *tau;
#line 320 "dlaed6.f"
    } else {
#line 321 "dlaed6.f"
	ubd = *tau;
#line 322 "dlaed6.f"
    }

/*        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent */
/*                            scheme */

/*     It is not hard to see that */

/*           1) Iterations will go up monotonically */
/*              if FINIT < 0; */

/*           2) Iterations will go down monotonically */
/*              if FINIT > 0. */

#line 335 "dlaed6.f"
    iter = niter + 1;

#line 337 "dlaed6.f"
    for (niter = iter; niter <= 40; ++niter) {

#line 339 "dlaed6.f"
	if (*orgati) {
#line 340 "dlaed6.f"
	    temp1 = dscale[1] - *tau;
#line 341 "dlaed6.f"
	    temp2 = dscale[2] - *tau;
#line 342 "dlaed6.f"
	} else {
#line 343 "dlaed6.f"
	    temp1 = dscale[0] - *tau;
#line 344 "dlaed6.f"
	    temp2 = dscale[1] - *tau;
#line 345 "dlaed6.f"
	}
#line 346 "dlaed6.f"
	a = (temp1 + temp2) * f - temp1 * temp2 * df;
#line 347 "dlaed6.f"
	b = temp1 * temp2 * f;
#line 348 "dlaed6.f"
	c__ = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
/* Computing MAX */
#line 349 "dlaed6.f"
	d__1 = abs(a), d__2 = abs(b), d__1 = max(d__1,d__2), d__2 = abs(c__);
#line 349 "dlaed6.f"
	temp = max(d__1,d__2);
#line 350 "dlaed6.f"
	a /= temp;
#line 351 "dlaed6.f"
	b /= temp;
#line 352 "dlaed6.f"
	c__ /= temp;
#line 353 "dlaed6.f"
	if (c__ == 0.) {
#line 354 "dlaed6.f"
	    eta = b / a;
#line 355 "dlaed6.f"
	} else if (a <= 0.) {
#line 356 "dlaed6.f"
	    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ 
		    * 2.);
#line 357 "dlaed6.f"
	} else {
#line 358 "dlaed6.f"
	    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))
		    );
#line 359 "dlaed6.f"
	}
#line 360 "dlaed6.f"
	if (f * eta >= 0.) {
#line 361 "dlaed6.f"
	    eta = -f / df;
#line 362 "dlaed6.f"
	}

#line 364 "dlaed6.f"
	*tau += eta;
#line 365 "dlaed6.f"
	if (*tau < lbd || *tau > ubd) {
#line 365 "dlaed6.f"
	    *tau = (lbd + ubd) / 2.;
#line 365 "dlaed6.f"
	}

#line 368 "dlaed6.f"
	fc = 0.;
#line 369 "dlaed6.f"
	erretm = 0.;
#line 370 "dlaed6.f"
	df = 0.;
#line 371 "dlaed6.f"
	ddf = 0.;
#line 372 "dlaed6.f"
	for (i__ = 1; i__ <= 3; ++i__) {
#line 373 "dlaed6.f"
	    if (dscale[i__ - 1] - *tau != 0.) {
#line 374 "dlaed6.f"
		temp = 1. / (dscale[i__ - 1] - *tau);
#line 375 "dlaed6.f"
		temp1 = zscale[i__ - 1] * temp;
#line 376 "dlaed6.f"
		temp2 = temp1 * temp;
#line 377 "dlaed6.f"
		temp3 = temp2 * temp;
#line 378 "dlaed6.f"
		temp4 = temp1 / dscale[i__ - 1];
#line 379 "dlaed6.f"
		fc += temp4;
#line 380 "dlaed6.f"
		erretm += abs(temp4);
#line 381 "dlaed6.f"
		df += temp2;
#line 382 "dlaed6.f"
		ddf += temp3;
#line 383 "dlaed6.f"
	    } else {
#line 384 "dlaed6.f"
		goto L60;
#line 385 "dlaed6.f"
	    }
#line 386 "dlaed6.f"
/* L40: */
#line 386 "dlaed6.f"
	}
#line 387 "dlaed6.f"
	f = *finit + *tau * fc;
#line 388 "dlaed6.f"
	erretm = (abs(*finit) + abs(*tau) * erretm) * 8. + abs(*tau) * df;
#line 390 "dlaed6.f"
	if (abs(f) <= eps * 4. * erretm || ubd - lbd <= eps * 4. * abs(*tau)) 
		{
#line 390 "dlaed6.f"
	    goto L60;
#line 390 "dlaed6.f"
	}
#line 393 "dlaed6.f"
	if (f <= 0.) {
#line 394 "dlaed6.f"
	    lbd = *tau;
#line 395 "dlaed6.f"
	} else {
#line 396 "dlaed6.f"
	    ubd = *tau;
#line 397 "dlaed6.f"
	}
#line 398 "dlaed6.f"
/* L50: */
#line 398 "dlaed6.f"
    }
#line 399 "dlaed6.f"
    *info = 1;
#line 400 "dlaed6.f"
L60:

/*     Undo scaling */

#line 404 "dlaed6.f"
    if (scale) {
#line 404 "dlaed6.f"
	*tau *= sclinv;
#line 404 "dlaed6.f"
    }
#line 406 "dlaed6.f"
    return 0;

/*     End of DLAED6 */

} /* dlaed6_ */

