#line 1 "dlasq4.f"
/* dlasq4.f -- translated by f2c (version 20100827).
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

#line 1 "dlasq4.f"
/* > \brief \b DLASQ4 computes an approximation to the smallest eigenvalue using values of d from the previous
 transform. Used by sbdsqr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASQ4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq4.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq4.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq4.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, */
/*                          DN1, DN2, TAU, TTYPE, G ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I0, N0, N0IN, PP, TTYPE */
/*       DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ4 computes an approximation TAU to the smallest eigenvalue */
/* > using values of d from the previous transform. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] I0 */
/* > \verbatim */
/* >          I0 is INTEGER */
/* >        First index. */
/* > \endverbatim */
/* > */
/* > \param[in] N0 */
/* > \verbatim */
/* >          N0 is INTEGER */
/* >        Last index. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( 4*N0 ) */
/* >        Z holds the qd array. */
/* > \endverbatim */
/* > */
/* > \param[in] PP */
/* > \verbatim */
/* >          PP is INTEGER */
/* >        PP=0 for ping, PP=1 for pong. */
/* > \endverbatim */
/* > */
/* > \param[in] N0IN */
/* > \verbatim */
/* >          N0IN is INTEGER */
/* >        The value of N0 at start of EIGTEST. */
/* > \endverbatim */
/* > */
/* > \param[in] DMIN */
/* > \verbatim */
/* >          DMIN is DOUBLE PRECISION */
/* >        Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[in] DMIN1 */
/* > \verbatim */
/* >          DMIN1 is DOUBLE PRECISION */
/* >        Minimum value of d, excluding D( N0 ). */
/* > \endverbatim */
/* > */
/* > \param[in] DMIN2 */
/* > \verbatim */
/* >          DMIN2 is DOUBLE PRECISION */
/* >        Minimum value of d, excluding D( N0 ) and D( N0-1 ). */
/* > \endverbatim */
/* > */
/* > \param[in] DN */
/* > \verbatim */
/* >          DN is DOUBLE PRECISION */
/* >        d(N) */
/* > \endverbatim */
/* > */
/* > \param[in] DN1 */
/* > \verbatim */
/* >          DN1 is DOUBLE PRECISION */
/* >        d(N-1) */
/* > \endverbatim */
/* > */
/* > \param[in] DN2 */
/* > \verbatim */
/* >          DN2 is DOUBLE PRECISION */
/* >        d(N-2) */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION */
/* >        This is the shift. */
/* > \endverbatim */
/* > */
/* > \param[out] TTYPE */
/* > \verbatim */
/* >          TTYPE is INTEGER */
/* >        Shift type. */
/* > \endverbatim */
/* > */
/* > \param[in,out] G */
/* > \verbatim */
/* >          G is DOUBLE PRECISION */
/* >        G is passed as an argument in order to save its value between */
/* >        calls to DLASQ4. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup auxOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  CNST1 = 9/16 */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlasq4_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, integer *n0in, doublereal *dmin__, doublereal *dmin1, 
	doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, 
	doublereal *tau, integer *ttype, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal s, a2, b1, b2;
    static integer i4, nn, np;
    static doublereal gam, gap1, gap2;


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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
/*     .. Executable Statements .. */

/*     A negative DMIN forces the shift to take that absolute value */
/*     TTYPE records the type of shift. */

#line 190 "dlasq4.f"
    /* Parameter adjustments */
#line 190 "dlasq4.f"
    --z__;
#line 190 "dlasq4.f"

#line 190 "dlasq4.f"
    /* Function Body */
#line 190 "dlasq4.f"
    if (*dmin__ <= 0.) {
#line 191 "dlasq4.f"
	*tau = -(*dmin__);
#line 192 "dlasq4.f"
	*ttype = -1;
#line 193 "dlasq4.f"
	return 0;
#line 194 "dlasq4.f"
    }

#line 196 "dlasq4.f"
    nn = (*n0 << 2) + *pp;
#line 197 "dlasq4.f"
    if (*n0in == *n0) {

/*        No eigenvalues deflated. */

#line 201 "dlasq4.f"
	if (*dmin__ == *dn || *dmin__ == *dn1) {

#line 203 "dlasq4.f"
	    b1 = sqrt(z__[nn - 3]) * sqrt(z__[nn - 5]);
#line 204 "dlasq4.f"
	    b2 = sqrt(z__[nn - 7]) * sqrt(z__[nn - 9]);
#line 205 "dlasq4.f"
	    a2 = z__[nn - 7] + z__[nn - 5];

/*           Cases 2 and 3. */

#line 209 "dlasq4.f"
	    if (*dmin__ == *dn && *dmin1 == *dn1) {
#line 210 "dlasq4.f"
		gap2 = *dmin2 - a2 - *dmin2 * .25;
#line 211 "dlasq4.f"
		if (gap2 > 0. && gap2 > b2) {
#line 212 "dlasq4.f"
		    gap1 = a2 - *dn - b2 / gap2 * b2;
#line 213 "dlasq4.f"
		} else {
#line 214 "dlasq4.f"
		    gap1 = a2 - *dn - (b1 + b2);
#line 215 "dlasq4.f"
		}
#line 216 "dlasq4.f"
		if (gap1 > 0. && gap1 > b1) {
/* Computing MAX */
#line 217 "dlasq4.f"
		    d__1 = *dn - b1 / gap1 * b1, d__2 = *dmin__ * .5;
#line 217 "dlasq4.f"
		    s = max(d__1,d__2);
#line 218 "dlasq4.f"
		    *ttype = -2;
#line 219 "dlasq4.f"
		} else {
#line 220 "dlasq4.f"
		    s = 0.;
#line 221 "dlasq4.f"
		    if (*dn > b1) {
#line 221 "dlasq4.f"
			s = *dn - b1;
#line 221 "dlasq4.f"
		    }
#line 223 "dlasq4.f"
		    if (a2 > b1 + b2) {
/* Computing MIN */
#line 223 "dlasq4.f"
			d__1 = s, d__2 = a2 - (b1 + b2);
#line 223 "dlasq4.f"
			s = min(d__1,d__2);
#line 223 "dlasq4.f"
		    }
/* Computing MAX */
#line 225 "dlasq4.f"
		    d__1 = s, d__2 = *dmin__ * .333;
#line 225 "dlasq4.f"
		    s = max(d__1,d__2);
#line 226 "dlasq4.f"
		    *ttype = -3;
#line 227 "dlasq4.f"
		}
#line 228 "dlasq4.f"
	    } else {

/*              Case 4. */

#line 232 "dlasq4.f"
		*ttype = -4;
#line 233 "dlasq4.f"
		s = *dmin__ * .25;
#line 234 "dlasq4.f"
		if (*dmin__ == *dn) {
#line 235 "dlasq4.f"
		    gam = *dn;
#line 236 "dlasq4.f"
		    a2 = 0.;
#line 237 "dlasq4.f"
		    if (z__[nn - 5] > z__[nn - 7]) {
#line 237 "dlasq4.f"
			return 0;
#line 237 "dlasq4.f"
		    }
#line 239 "dlasq4.f"
		    b2 = z__[nn - 5] / z__[nn - 7];
#line 240 "dlasq4.f"
		    np = nn - 9;
#line 241 "dlasq4.f"
		} else {
#line 242 "dlasq4.f"
		    np = nn - (*pp << 1);
#line 243 "dlasq4.f"
		    gam = *dn1;
#line 244 "dlasq4.f"
		    if (z__[np - 4] > z__[np - 2]) {
#line 244 "dlasq4.f"
			return 0;
#line 244 "dlasq4.f"
		    }
#line 246 "dlasq4.f"
		    a2 = z__[np - 4] / z__[np - 2];
#line 247 "dlasq4.f"
		    if (z__[nn - 9] > z__[nn - 11]) {
#line 247 "dlasq4.f"
			return 0;
#line 247 "dlasq4.f"
		    }
#line 249 "dlasq4.f"
		    b2 = z__[nn - 9] / z__[nn - 11];
#line 250 "dlasq4.f"
		    np = nn - 13;
#line 251 "dlasq4.f"
		}

/*              Approximate contribution to norm squared from I < NN-1. */

#line 255 "dlasq4.f"
		a2 += b2;
#line 256 "dlasq4.f"
		i__1 = (*i0 << 2) - 1 + *pp;
#line 256 "dlasq4.f"
		for (i4 = np; i4 >= i__1; i4 += -4) {
#line 257 "dlasq4.f"
		    if (b2 == 0.) {
#line 257 "dlasq4.f"
			goto L20;
#line 257 "dlasq4.f"
		    }
#line 259 "dlasq4.f"
		    b1 = b2;
#line 260 "dlasq4.f"
		    if (z__[i4] > z__[i4 - 2]) {
#line 260 "dlasq4.f"
			return 0;
#line 260 "dlasq4.f"
		    }
#line 262 "dlasq4.f"
		    b2 *= z__[i4] / z__[i4 - 2];
#line 263 "dlasq4.f"
		    a2 += b2;
#line 264 "dlasq4.f"
		    if (max(b2,b1) * 100. < a2 || .563 < a2) {
#line 264 "dlasq4.f"
			goto L20;
#line 264 "dlasq4.f"
		    }
#line 266 "dlasq4.f"
/* L10: */
#line 266 "dlasq4.f"
		}
#line 267 "dlasq4.f"
L20:
#line 268 "dlasq4.f"
		a2 *= 1.05;

/*              Rayleigh quotient residual bound. */

#line 272 "dlasq4.f"
		if (a2 < .563) {
#line 272 "dlasq4.f"
		    s = gam * (1. - sqrt(a2)) / (a2 + 1.);
#line 272 "dlasq4.f"
		}
#line 274 "dlasq4.f"
	    }
#line 275 "dlasq4.f"
	} else if (*dmin__ == *dn2) {

/*           Case 5. */

#line 279 "dlasq4.f"
	    *ttype = -5;
#line 280 "dlasq4.f"
	    s = *dmin__ * .25;

/*           Compute contribution to norm squared from I > NN-2. */

#line 284 "dlasq4.f"
	    np = nn - (*pp << 1);
#line 285 "dlasq4.f"
	    b1 = z__[np - 2];
#line 286 "dlasq4.f"
	    b2 = z__[np - 6];
#line 287 "dlasq4.f"
	    gam = *dn2;
#line 288 "dlasq4.f"
	    if (z__[np - 8] > b2 || z__[np - 4] > b1) {
#line 288 "dlasq4.f"
		return 0;
#line 288 "dlasq4.f"
	    }
#line 290 "dlasq4.f"
	    a2 = z__[np - 8] / b2 * (z__[np - 4] / b1 + 1.);

/*           Approximate contribution to norm squared from I < NN-2. */

#line 294 "dlasq4.f"
	    if (*n0 - *i0 > 2) {
#line 295 "dlasq4.f"
		b2 = z__[nn - 13] / z__[nn - 15];
#line 296 "dlasq4.f"
		a2 += b2;
#line 297 "dlasq4.f"
		i__1 = (*i0 << 2) - 1 + *pp;
#line 297 "dlasq4.f"
		for (i4 = nn - 17; i4 >= i__1; i4 += -4) {
#line 298 "dlasq4.f"
		    if (b2 == 0.) {
#line 298 "dlasq4.f"
			goto L40;
#line 298 "dlasq4.f"
		    }
#line 300 "dlasq4.f"
		    b1 = b2;
#line 301 "dlasq4.f"
		    if (z__[i4] > z__[i4 - 2]) {
#line 301 "dlasq4.f"
			return 0;
#line 301 "dlasq4.f"
		    }
#line 303 "dlasq4.f"
		    b2 *= z__[i4] / z__[i4 - 2];
#line 304 "dlasq4.f"
		    a2 += b2;
#line 305 "dlasq4.f"
		    if (max(b2,b1) * 100. < a2 || .563 < a2) {
#line 305 "dlasq4.f"
			goto L40;
#line 305 "dlasq4.f"
		    }
#line 307 "dlasq4.f"
/* L30: */
#line 307 "dlasq4.f"
		}
#line 308 "dlasq4.f"
L40:
#line 309 "dlasq4.f"
		a2 *= 1.05;
#line 310 "dlasq4.f"
	    }

#line 312 "dlasq4.f"
	    if (a2 < .563) {
#line 312 "dlasq4.f"
		s = gam * (1. - sqrt(a2)) / (a2 + 1.);
#line 312 "dlasq4.f"
	    }
#line 314 "dlasq4.f"
	} else {

/*           Case 6, no information to guide us. */

#line 318 "dlasq4.f"
	    if (*ttype == -6) {
#line 319 "dlasq4.f"
		*g += (1. - *g) * .333;
#line 320 "dlasq4.f"
	    } else if (*ttype == -18) {
#line 321 "dlasq4.f"
		*g = .083250000000000005;
#line 322 "dlasq4.f"
	    } else {
#line 323 "dlasq4.f"
		*g = .25;
#line 324 "dlasq4.f"
	    }
#line 325 "dlasq4.f"
	    s = *g * *dmin__;
#line 326 "dlasq4.f"
	    *ttype = -6;
#line 327 "dlasq4.f"
	}

#line 329 "dlasq4.f"
    } else if (*n0in == *n0 + 1) {

/*        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN. */

#line 333 "dlasq4.f"
	if (*dmin1 == *dn1 && *dmin2 == *dn2) {

/*           Cases 7 and 8. */

#line 337 "dlasq4.f"
	    *ttype = -7;
#line 338 "dlasq4.f"
	    s = *dmin1 * .333;
#line 339 "dlasq4.f"
	    if (z__[nn - 5] > z__[nn - 7]) {
#line 339 "dlasq4.f"
		return 0;
#line 339 "dlasq4.f"
	    }
#line 341 "dlasq4.f"
	    b1 = z__[nn - 5] / z__[nn - 7];
#line 342 "dlasq4.f"
	    b2 = b1;
#line 343 "dlasq4.f"
	    if (b2 == 0.) {
#line 343 "dlasq4.f"
		goto L60;
#line 343 "dlasq4.f"
	    }
#line 345 "dlasq4.f"
	    i__1 = (*i0 << 2) - 1 + *pp;
#line 345 "dlasq4.f"
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
#line 346 "dlasq4.f"
		a2 = b1;
#line 347 "dlasq4.f"
		if (z__[i4] > z__[i4 - 2]) {
#line 347 "dlasq4.f"
		    return 0;
#line 347 "dlasq4.f"
		}
#line 349 "dlasq4.f"
		b1 *= z__[i4] / z__[i4 - 2];
#line 350 "dlasq4.f"
		b2 += b1;
#line 351 "dlasq4.f"
		if (max(b1,a2) * 100. < b2) {
#line 351 "dlasq4.f"
		    goto L60;
#line 351 "dlasq4.f"
		}
#line 353 "dlasq4.f"
/* L50: */
#line 353 "dlasq4.f"
	    }
#line 354 "dlasq4.f"
L60:
#line 355 "dlasq4.f"
	    b2 = sqrt(b2 * 1.05);
/* Computing 2nd power */
#line 356 "dlasq4.f"
	    d__1 = b2;
#line 356 "dlasq4.f"
	    a2 = *dmin1 / (d__1 * d__1 + 1.);
#line 357 "dlasq4.f"
	    gap2 = *dmin2 * .5 - a2;
#line 358 "dlasq4.f"
	    if (gap2 > 0. && gap2 > b2 * a2) {
/* Computing MAX */
#line 359 "dlasq4.f"
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
#line 359 "dlasq4.f"
		s = max(d__1,d__2);
#line 360 "dlasq4.f"
	    } else {
/* Computing MAX */
#line 361 "dlasq4.f"
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
#line 361 "dlasq4.f"
		s = max(d__1,d__2);
#line 362 "dlasq4.f"
		*ttype = -8;
#line 363 "dlasq4.f"
	    }
#line 364 "dlasq4.f"
	} else {

/*           Case 9. */

#line 368 "dlasq4.f"
	    s = *dmin1 * .25;
#line 369 "dlasq4.f"
	    if (*dmin1 == *dn1) {
#line 369 "dlasq4.f"
		s = *dmin1 * .5;
#line 369 "dlasq4.f"
	    }
#line 371 "dlasq4.f"
	    *ttype = -9;
#line 372 "dlasq4.f"
	}

#line 374 "dlasq4.f"
    } else if (*n0in == *n0 + 2) {

/*        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN. */

/*        Cases 10 and 11. */

#line 380 "dlasq4.f"
	if (*dmin2 == *dn2 && z__[nn - 5] * 2. < z__[nn - 7]) {
#line 381 "dlasq4.f"
	    *ttype = -10;
#line 382 "dlasq4.f"
	    s = *dmin2 * .333;
#line 383 "dlasq4.f"
	    if (z__[nn - 5] > z__[nn - 7]) {
#line 383 "dlasq4.f"
		return 0;
#line 383 "dlasq4.f"
	    }
#line 385 "dlasq4.f"
	    b1 = z__[nn - 5] / z__[nn - 7];
#line 386 "dlasq4.f"
	    b2 = b1;
#line 387 "dlasq4.f"
	    if (b2 == 0.) {
#line 387 "dlasq4.f"
		goto L80;
#line 387 "dlasq4.f"
	    }
#line 389 "dlasq4.f"
	    i__1 = (*i0 << 2) - 1 + *pp;
#line 389 "dlasq4.f"
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
#line 390 "dlasq4.f"
		if (z__[i4] > z__[i4 - 2]) {
#line 390 "dlasq4.f"
		    return 0;
#line 390 "dlasq4.f"
		}
#line 392 "dlasq4.f"
		b1 *= z__[i4] / z__[i4 - 2];
#line 393 "dlasq4.f"
		b2 += b1;
#line 394 "dlasq4.f"
		if (b1 * 100. < b2) {
#line 394 "dlasq4.f"
		    goto L80;
#line 394 "dlasq4.f"
		}
#line 396 "dlasq4.f"
/* L70: */
#line 396 "dlasq4.f"
	    }
#line 397 "dlasq4.f"
L80:
#line 398 "dlasq4.f"
	    b2 = sqrt(b2 * 1.05);
/* Computing 2nd power */
#line 399 "dlasq4.f"
	    d__1 = b2;
#line 399 "dlasq4.f"
	    a2 = *dmin2 / (d__1 * d__1 + 1.);
#line 400 "dlasq4.f"
	    gap2 = z__[nn - 7] + z__[nn - 9] - sqrt(z__[nn - 11]) * sqrt(z__[
		    nn - 9]) - a2;
#line 402 "dlasq4.f"
	    if (gap2 > 0. && gap2 > b2 * a2) {
/* Computing MAX */
#line 403 "dlasq4.f"
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
#line 403 "dlasq4.f"
		s = max(d__1,d__2);
#line 404 "dlasq4.f"
	    } else {
/* Computing MAX */
#line 405 "dlasq4.f"
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
#line 405 "dlasq4.f"
		s = max(d__1,d__2);
#line 406 "dlasq4.f"
	    }
#line 407 "dlasq4.f"
	} else {
#line 408 "dlasq4.f"
	    s = *dmin2 * .25;
#line 409 "dlasq4.f"
	    *ttype = -11;
#line 410 "dlasq4.f"
	}
#line 411 "dlasq4.f"
    } else if (*n0in > *n0 + 2) {

/*        Case 12, more than two eigenvalues deflated. No information. */

#line 415 "dlasq4.f"
	s = 0.;
#line 416 "dlasq4.f"
	*ttype = -12;
#line 417 "dlasq4.f"
    }

#line 419 "dlasq4.f"
    *tau = s;
#line 420 "dlasq4.f"
    return 0;

/*     End of DLASQ4 */

} /* dlasq4_ */

