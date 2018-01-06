#line 1 "slasq4.f"
/* slasq4.f -- translated by f2c (version 20100827).
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

#line 1 "slasq4.f"
/* > \brief \b SLASQ4 computes an approximation to the smallest eigenvalue using values of d from the previous
 transform. Used by sbdsqr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASQ4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq4.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq4.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq4.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, */
/*                          DN1, DN2, TAU, TTYPE, G ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I0, N0, N0IN, PP, TTYPE */
/*       REAL               DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASQ4 computes an approximation TAU to the smallest eigenvalue */
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
/* >          Z is REAL array, dimension ( 4*N ) */
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
/* >          DMIN is REAL */
/* >        Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[in] DMIN1 */
/* > \verbatim */
/* >          DMIN1 is REAL */
/* >        Minimum value of d, excluding D( N0 ). */
/* > \endverbatim */
/* > */
/* > \param[in] DMIN2 */
/* > \verbatim */
/* >          DMIN2 is REAL */
/* >        Minimum value of d, excluding D( N0 ) and D( N0-1 ). */
/* > \endverbatim */
/* > */
/* > \param[in] DN */
/* > \verbatim */
/* >          DN is REAL */
/* >        d(N) */
/* > \endverbatim */
/* > */
/* > \param[in] DN1 */
/* > \verbatim */
/* >          DN1 is REAL */
/* >        d(N-1) */
/* > \endverbatim */
/* > */
/* > \param[in] DN2 */
/* > \verbatim */
/* >          DN2 is REAL */
/* >        d(N-2) */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL */
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
/* >          G is REAL */
/* >        G is passed as an argument in order to save its value between */
/* >        calls to SLASQ4. */
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
/* >  CNST1 = 9/16 */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slasq4_(integer *i0, integer *n0, doublereal *z__, 
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
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     A negative DMIN forces the shift to take that absolute value */
/*     TTYPE records the type of shift. */

#line 190 "slasq4.f"
    /* Parameter adjustments */
#line 190 "slasq4.f"
    --z__;
#line 190 "slasq4.f"

#line 190 "slasq4.f"
    /* Function Body */
#line 190 "slasq4.f"
    if (*dmin__ <= 0.) {
#line 191 "slasq4.f"
	*tau = -(*dmin__);
#line 192 "slasq4.f"
	*ttype = -1;
#line 193 "slasq4.f"
	return 0;
#line 194 "slasq4.f"
    }

#line 196 "slasq4.f"
    nn = (*n0 << 2) + *pp;
#line 197 "slasq4.f"
    if (*n0in == *n0) {

/*        No eigenvalues deflated. */

#line 201 "slasq4.f"
	if (*dmin__ == *dn || *dmin__ == *dn1) {

#line 203 "slasq4.f"
	    b1 = sqrt(z__[nn - 3]) * sqrt(z__[nn - 5]);
#line 204 "slasq4.f"
	    b2 = sqrt(z__[nn - 7]) * sqrt(z__[nn - 9]);
#line 205 "slasq4.f"
	    a2 = z__[nn - 7] + z__[nn - 5];

/*           Cases 2 and 3. */

#line 209 "slasq4.f"
	    if (*dmin__ == *dn && *dmin1 == *dn1) {
#line 210 "slasq4.f"
		gap2 = *dmin2 - a2 - *dmin2 * .25;
#line 211 "slasq4.f"
		if (gap2 > 0. && gap2 > b2) {
#line 212 "slasq4.f"
		    gap1 = a2 - *dn - b2 / gap2 * b2;
#line 213 "slasq4.f"
		} else {
#line 214 "slasq4.f"
		    gap1 = a2 - *dn - (b1 + b2);
#line 215 "slasq4.f"
		}
#line 216 "slasq4.f"
		if (gap1 > 0. && gap1 > b1) {
/* Computing MAX */
#line 217 "slasq4.f"
		    d__1 = *dn - b1 / gap1 * b1, d__2 = *dmin__ * .5;
#line 217 "slasq4.f"
		    s = max(d__1,d__2);
#line 218 "slasq4.f"
		    *ttype = -2;
#line 219 "slasq4.f"
		} else {
#line 220 "slasq4.f"
		    s = 0.;
#line 221 "slasq4.f"
		    if (*dn > b1) {
#line 221 "slasq4.f"
			s = *dn - b1;
#line 221 "slasq4.f"
		    }
#line 223 "slasq4.f"
		    if (a2 > b1 + b2) {
/* Computing MIN */
#line 223 "slasq4.f"
			d__1 = s, d__2 = a2 - (b1 + b2);
#line 223 "slasq4.f"
			s = min(d__1,d__2);
#line 223 "slasq4.f"
		    }
/* Computing MAX */
#line 225 "slasq4.f"
		    d__1 = s, d__2 = *dmin__ * .333;
#line 225 "slasq4.f"
		    s = max(d__1,d__2);
#line 226 "slasq4.f"
		    *ttype = -3;
#line 227 "slasq4.f"
		}
#line 228 "slasq4.f"
	    } else {

/*              Case 4. */

#line 232 "slasq4.f"
		*ttype = -4;
#line 233 "slasq4.f"
		s = *dmin__ * .25;
#line 234 "slasq4.f"
		if (*dmin__ == *dn) {
#line 235 "slasq4.f"
		    gam = *dn;
#line 236 "slasq4.f"
		    a2 = 0.;
#line 237 "slasq4.f"
		    if (z__[nn - 5] > z__[nn - 7]) {
#line 237 "slasq4.f"
			return 0;
#line 237 "slasq4.f"
		    }
#line 239 "slasq4.f"
		    b2 = z__[nn - 5] / z__[nn - 7];
#line 240 "slasq4.f"
		    np = nn - 9;
#line 241 "slasq4.f"
		} else {
#line 242 "slasq4.f"
		    np = nn - (*pp << 1);
#line 243 "slasq4.f"
		    b2 = z__[np - 2];
#line 244 "slasq4.f"
		    gam = *dn1;
#line 245 "slasq4.f"
		    if (z__[np - 4] > z__[np - 2]) {
#line 245 "slasq4.f"
			return 0;
#line 245 "slasq4.f"
		    }
#line 247 "slasq4.f"
		    a2 = z__[np - 4] / z__[np - 2];
#line 248 "slasq4.f"
		    if (z__[nn - 9] > z__[nn - 11]) {
#line 248 "slasq4.f"
			return 0;
#line 248 "slasq4.f"
		    }
#line 250 "slasq4.f"
		    b2 = z__[nn - 9] / z__[nn - 11];
#line 251 "slasq4.f"
		    np = nn - 13;
#line 252 "slasq4.f"
		}

/*              Approximate contribution to norm squared from I < NN-1. */

#line 256 "slasq4.f"
		a2 += b2;
#line 257 "slasq4.f"
		i__1 = (*i0 << 2) - 1 + *pp;
#line 257 "slasq4.f"
		for (i4 = np; i4 >= i__1; i4 += -4) {
#line 258 "slasq4.f"
		    if (b2 == 0.) {
#line 258 "slasq4.f"
			goto L20;
#line 258 "slasq4.f"
		    }
#line 260 "slasq4.f"
		    b1 = b2;
#line 261 "slasq4.f"
		    if (z__[i4] > z__[i4 - 2]) {
#line 261 "slasq4.f"
			return 0;
#line 261 "slasq4.f"
		    }
#line 263 "slasq4.f"
		    b2 *= z__[i4] / z__[i4 - 2];
#line 264 "slasq4.f"
		    a2 += b2;
#line 265 "slasq4.f"
		    if (max(b2,b1) * 100. < a2 || .563 < a2) {
#line 265 "slasq4.f"
			goto L20;
#line 265 "slasq4.f"
		    }
#line 267 "slasq4.f"
/* L10: */
#line 267 "slasq4.f"
		}
#line 268 "slasq4.f"
L20:
#line 269 "slasq4.f"
		a2 *= 1.05;

/*              Rayleigh quotient residual bound. */

#line 273 "slasq4.f"
		if (a2 < .563) {
#line 273 "slasq4.f"
		    s = gam * (1. - sqrt(a2)) / (a2 + 1.);
#line 273 "slasq4.f"
		}
#line 275 "slasq4.f"
	    }
#line 276 "slasq4.f"
	} else if (*dmin__ == *dn2) {

/*           Case 5. */

#line 280 "slasq4.f"
	    *ttype = -5;
#line 281 "slasq4.f"
	    s = *dmin__ * .25;

/*           Compute contribution to norm squared from I > NN-2. */

#line 285 "slasq4.f"
	    np = nn - (*pp << 1);
#line 286 "slasq4.f"
	    b1 = z__[np - 2];
#line 287 "slasq4.f"
	    b2 = z__[np - 6];
#line 288 "slasq4.f"
	    gam = *dn2;
#line 289 "slasq4.f"
	    if (z__[np - 8] > b2 || z__[np - 4] > b1) {
#line 289 "slasq4.f"
		return 0;
#line 289 "slasq4.f"
	    }
#line 291 "slasq4.f"
	    a2 = z__[np - 8] / b2 * (z__[np - 4] / b1 + 1.);

/*           Approximate contribution to norm squared from I < NN-2. */

#line 295 "slasq4.f"
	    if (*n0 - *i0 > 2) {
#line 296 "slasq4.f"
		b2 = z__[nn - 13] / z__[nn - 15];
#line 297 "slasq4.f"
		a2 += b2;
#line 298 "slasq4.f"
		i__1 = (*i0 << 2) - 1 + *pp;
#line 298 "slasq4.f"
		for (i4 = nn - 17; i4 >= i__1; i4 += -4) {
#line 299 "slasq4.f"
		    if (b2 == 0.) {
#line 299 "slasq4.f"
			goto L40;
#line 299 "slasq4.f"
		    }
#line 301 "slasq4.f"
		    b1 = b2;
#line 302 "slasq4.f"
		    if (z__[i4] > z__[i4 - 2]) {
#line 302 "slasq4.f"
			return 0;
#line 302 "slasq4.f"
		    }
#line 304 "slasq4.f"
		    b2 *= z__[i4] / z__[i4 - 2];
#line 305 "slasq4.f"
		    a2 += b2;
#line 306 "slasq4.f"
		    if (max(b2,b1) * 100. < a2 || .563 < a2) {
#line 306 "slasq4.f"
			goto L40;
#line 306 "slasq4.f"
		    }
#line 308 "slasq4.f"
/* L30: */
#line 308 "slasq4.f"
		}
#line 309 "slasq4.f"
L40:
#line 310 "slasq4.f"
		a2 *= 1.05;
#line 311 "slasq4.f"
	    }

#line 313 "slasq4.f"
	    if (a2 < .563) {
#line 313 "slasq4.f"
		s = gam * (1. - sqrt(a2)) / (a2 + 1.);
#line 313 "slasq4.f"
	    }
#line 315 "slasq4.f"
	} else {

/*           Case 6, no information to guide us. */

#line 319 "slasq4.f"
	    if (*ttype == -6) {
#line 320 "slasq4.f"
		*g += (1. - *g) * .333;
#line 321 "slasq4.f"
	    } else if (*ttype == -18) {
#line 322 "slasq4.f"
		*g = .083250000000000005;
#line 323 "slasq4.f"
	    } else {
#line 324 "slasq4.f"
		*g = .25;
#line 325 "slasq4.f"
	    }
#line 326 "slasq4.f"
	    s = *g * *dmin__;
#line 327 "slasq4.f"
	    *ttype = -6;
#line 328 "slasq4.f"
	}

#line 330 "slasq4.f"
    } else if (*n0in == *n0 + 1) {

/*        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN. */

#line 334 "slasq4.f"
	if (*dmin1 == *dn1 && *dmin2 == *dn2) {

/*           Cases 7 and 8. */

#line 338 "slasq4.f"
	    *ttype = -7;
#line 339 "slasq4.f"
	    s = *dmin1 * .333;
#line 340 "slasq4.f"
	    if (z__[nn - 5] > z__[nn - 7]) {
#line 340 "slasq4.f"
		return 0;
#line 340 "slasq4.f"
	    }
#line 342 "slasq4.f"
	    b1 = z__[nn - 5] / z__[nn - 7];
#line 343 "slasq4.f"
	    b2 = b1;
#line 344 "slasq4.f"
	    if (b2 == 0.) {
#line 344 "slasq4.f"
		goto L60;
#line 344 "slasq4.f"
	    }
#line 346 "slasq4.f"
	    i__1 = (*i0 << 2) - 1 + *pp;
#line 346 "slasq4.f"
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
#line 347 "slasq4.f"
		a2 = b1;
#line 348 "slasq4.f"
		if (z__[i4] > z__[i4 - 2]) {
#line 348 "slasq4.f"
		    return 0;
#line 348 "slasq4.f"
		}
#line 350 "slasq4.f"
		b1 *= z__[i4] / z__[i4 - 2];
#line 351 "slasq4.f"
		b2 += b1;
#line 352 "slasq4.f"
		if (max(b1,a2) * 100. < b2) {
#line 352 "slasq4.f"
		    goto L60;
#line 352 "slasq4.f"
		}
#line 354 "slasq4.f"
/* L50: */
#line 354 "slasq4.f"
	    }
#line 355 "slasq4.f"
L60:
#line 356 "slasq4.f"
	    b2 = sqrt(b2 * 1.05);
/* Computing 2nd power */
#line 357 "slasq4.f"
	    d__1 = b2;
#line 357 "slasq4.f"
	    a2 = *dmin1 / (d__1 * d__1 + 1.);
#line 358 "slasq4.f"
	    gap2 = *dmin2 * .5 - a2;
#line 359 "slasq4.f"
	    if (gap2 > 0. && gap2 > b2 * a2) {
/* Computing MAX */
#line 360 "slasq4.f"
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
#line 360 "slasq4.f"
		s = max(d__1,d__2);
#line 361 "slasq4.f"
	    } else {
/* Computing MAX */
#line 362 "slasq4.f"
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
#line 362 "slasq4.f"
		s = max(d__1,d__2);
#line 363 "slasq4.f"
		*ttype = -8;
#line 364 "slasq4.f"
	    }
#line 365 "slasq4.f"
	} else {

/*           Case 9. */

#line 369 "slasq4.f"
	    s = *dmin1 * .25;
#line 370 "slasq4.f"
	    if (*dmin1 == *dn1) {
#line 370 "slasq4.f"
		s = *dmin1 * .5;
#line 370 "slasq4.f"
	    }
#line 372 "slasq4.f"
	    *ttype = -9;
#line 373 "slasq4.f"
	}

#line 375 "slasq4.f"
    } else if (*n0in == *n0 + 2) {

/*        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN. */

/*        Cases 10 and 11. */

#line 381 "slasq4.f"
	if (*dmin2 == *dn2 && z__[nn - 5] * 2. < z__[nn - 7]) {
#line 382 "slasq4.f"
	    *ttype = -10;
#line 383 "slasq4.f"
	    s = *dmin2 * .333;
#line 384 "slasq4.f"
	    if (z__[nn - 5] > z__[nn - 7]) {
#line 384 "slasq4.f"
		return 0;
#line 384 "slasq4.f"
	    }
#line 386 "slasq4.f"
	    b1 = z__[nn - 5] / z__[nn - 7];
#line 387 "slasq4.f"
	    b2 = b1;
#line 388 "slasq4.f"
	    if (b2 == 0.) {
#line 388 "slasq4.f"
		goto L80;
#line 388 "slasq4.f"
	    }
#line 390 "slasq4.f"
	    i__1 = (*i0 << 2) - 1 + *pp;
#line 390 "slasq4.f"
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
#line 391 "slasq4.f"
		if (z__[i4] > z__[i4 - 2]) {
#line 391 "slasq4.f"
		    return 0;
#line 391 "slasq4.f"
		}
#line 393 "slasq4.f"
		b1 *= z__[i4] / z__[i4 - 2];
#line 394 "slasq4.f"
		b2 += b1;
#line 395 "slasq4.f"
		if (b1 * 100. < b2) {
#line 395 "slasq4.f"
		    goto L80;
#line 395 "slasq4.f"
		}
#line 397 "slasq4.f"
/* L70: */
#line 397 "slasq4.f"
	    }
#line 398 "slasq4.f"
L80:
#line 399 "slasq4.f"
	    b2 = sqrt(b2 * 1.05);
/* Computing 2nd power */
#line 400 "slasq4.f"
	    d__1 = b2;
#line 400 "slasq4.f"
	    a2 = *dmin2 / (d__1 * d__1 + 1.);
#line 401 "slasq4.f"
	    gap2 = z__[nn - 7] + z__[nn - 9] - sqrt(z__[nn - 11]) * sqrt(z__[
		    nn - 9]) - a2;
#line 403 "slasq4.f"
	    if (gap2 > 0. && gap2 > b2 * a2) {
/* Computing MAX */
#line 404 "slasq4.f"
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
#line 404 "slasq4.f"
		s = max(d__1,d__2);
#line 405 "slasq4.f"
	    } else {
/* Computing MAX */
#line 406 "slasq4.f"
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
#line 406 "slasq4.f"
		s = max(d__1,d__2);
#line 407 "slasq4.f"
	    }
#line 408 "slasq4.f"
	} else {
#line 409 "slasq4.f"
	    s = *dmin2 * .25;
#line 410 "slasq4.f"
	    *ttype = -11;
#line 411 "slasq4.f"
	}
#line 412 "slasq4.f"
    } else if (*n0in > *n0 + 2) {

/*        Case 12, more than two eigenvalues deflated. No information. */

#line 416 "slasq4.f"
	s = 0.;
#line 417 "slasq4.f"
	*ttype = -12;
#line 418 "slasq4.f"
    }

#line 420 "slasq4.f"
    *tau = s;
#line 421 "slasq4.f"
    return 0;

/*     End of SLASQ4 */

} /* slasq4_ */

