#line 1 "slasq5.f"
/* slasq5.f -- translated by f2c (version 20100827).
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

#line 1 "slasq5.f"
/* > \brief \b SLASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASQ5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq5.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq5.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq5.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN, */
/*                          DNM1, DNM2, IEEE ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            IEEE */
/*       INTEGER            I0, N0, PP */
/*       REAL               DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASQ5 computes one dqds transform in ping-pong form, one */
/* > version for IEEE machines another for non IEEE machines. */
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
/* >        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid */
/* >        an extra argument. */
/* > \endverbatim */
/* > */
/* > \param[in] PP */
/* > \verbatim */
/* >          PP is INTEGER */
/* >        PP=0 for ping, PP=1 for pong. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is REAL */
/* >        This is the shift. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN */
/* > \verbatim */
/* >          DMIN is REAL */
/* >        Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN1 */
/* > \verbatim */
/* >          DMIN1 is REAL */
/* >        Minimum value of d, excluding D( N0 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN2 */
/* > \verbatim */
/* >          DMIN2 is REAL */
/* >        Minimum value of d, excluding D( N0 ) and D( N0-1 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DN */
/* > \verbatim */
/* >          DN is REAL */
/* >        d(N0), the last value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DNM1 */
/* > \verbatim */
/* >          DNM1 is REAL */
/* >        d(N0-1). */
/* > \endverbatim */
/* > */
/* > \param[out] DNM2 */
/* > \verbatim */
/* >          DNM2 is REAL */
/* >        d(N0-2). */
/* > \endverbatim */
/* > */
/* > \param[in] IEEE */
/* > \verbatim */
/* >          IEEE is LOGICAL */
/* >        Flag for IEEE or non IEEE arithmetic. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int slasq5_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *tau, doublereal *sigma, doublereal *dmin__, 
	doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *
	dnm1, doublereal *dnm2, logical *ieee, doublereal *eps)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer j4, j4p2;
    static doublereal emin, temp, dthresh;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameter .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 165 "slasq5.f"
    /* Parameter adjustments */
#line 165 "slasq5.f"
    --z__;
#line 165 "slasq5.f"

#line 165 "slasq5.f"
    /* Function Body */
#line 165 "slasq5.f"
    if (*n0 - *i0 - 1 <= 0) {
#line 165 "slasq5.f"
	return 0;
#line 165 "slasq5.f"
    }

#line 168 "slasq5.f"
    dthresh = *eps * (*sigma + *tau);
#line 169 "slasq5.f"
    if (*tau < dthresh * .5) {
#line 169 "slasq5.f"
	*tau = 0.;
#line 169 "slasq5.f"
    }
#line 170 "slasq5.f"
    if (*tau != 0.) {
#line 171 "slasq5.f"
	j4 = (*i0 << 2) + *pp - 3;
#line 172 "slasq5.f"
	emin = z__[j4 + 4];
#line 173 "slasq5.f"
	d__ = z__[j4] - *tau;
#line 174 "slasq5.f"
	*dmin__ = d__;
#line 175 "slasq5.f"
	*dmin1 = -z__[j4];

#line 177 "slasq5.f"
	if (*ieee) {

/*     Code for IEEE arithmetic. */

#line 181 "slasq5.f"
	    if (*pp == 0) {
#line 182 "slasq5.f"
		i__1 = *n0 - 3 << 2;
#line 182 "slasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 183 "slasq5.f"
		    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 184 "slasq5.f"
		    temp = z__[j4 + 1] / z__[j4 - 2];
#line 185 "slasq5.f"
		    d__ = d__ * temp - *tau;
#line 186 "slasq5.f"
		    *dmin__ = min(*dmin__,d__);
#line 187 "slasq5.f"
		    z__[j4] = z__[j4 - 1] * temp;
/* Computing MIN */
#line 188 "slasq5.f"
		    d__1 = z__[j4];
#line 188 "slasq5.f"
		    emin = min(d__1,emin);
#line 189 "slasq5.f"
/* L10: */
#line 189 "slasq5.f"
		}
#line 190 "slasq5.f"
	    } else {
#line 191 "slasq5.f"
		i__1 = *n0 - 3 << 2;
#line 191 "slasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 192 "slasq5.f"
		    z__[j4 - 3] = d__ + z__[j4];
#line 193 "slasq5.f"
		    temp = z__[j4 + 2] / z__[j4 - 3];
#line 194 "slasq5.f"
		    d__ = d__ * temp - *tau;
#line 195 "slasq5.f"
		    *dmin__ = min(*dmin__,d__);
#line 196 "slasq5.f"
		    z__[j4 - 1] = z__[j4] * temp;
/* Computing MIN */
#line 197 "slasq5.f"
		    d__1 = z__[j4 - 1];
#line 197 "slasq5.f"
		    emin = min(d__1,emin);
#line 198 "slasq5.f"
/* L20: */
#line 198 "slasq5.f"
		}
#line 199 "slasq5.f"
	    }

/*     Unroll last two steps. */

#line 203 "slasq5.f"
	    *dnm2 = d__;
#line 204 "slasq5.f"
	    *dmin2 = *dmin__;
#line 205 "slasq5.f"
	    j4 = (*n0 - 2 << 2) - *pp;
#line 206 "slasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 207 "slasq5.f"
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 208 "slasq5.f"
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 209 "slasq5.f"
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
#line 210 "slasq5.f"
	    *dmin__ = min(*dmin__,*dnm1);

#line 212 "slasq5.f"
	    *dmin1 = *dmin__;
#line 213 "slasq5.f"
	    j4 += 4;
#line 214 "slasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 215 "slasq5.f"
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 216 "slasq5.f"
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 217 "slasq5.f"
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
#line 218 "slasq5.f"
	    *dmin__ = min(*dmin__,*dn);

#line 220 "slasq5.f"
	} else {

/*     Code for non IEEE arithmetic. */

#line 224 "slasq5.f"
	    if (*pp == 0) {
#line 225 "slasq5.f"
		i__1 = *n0 - 3 << 2;
#line 225 "slasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 226 "slasq5.f"
		    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 227 "slasq5.f"
		    if (d__ < 0.) {
#line 228 "slasq5.f"
			return 0;
#line 229 "slasq5.f"
		    } else {
#line 230 "slasq5.f"
			z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
#line 231 "slasq5.f"
			d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
#line 232 "slasq5.f"
		    }
#line 233 "slasq5.f"
		    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 234 "slasq5.f"
		    d__1 = emin, d__2 = z__[j4];
#line 234 "slasq5.f"
		    emin = min(d__1,d__2);
#line 235 "slasq5.f"
/* L30: */
#line 235 "slasq5.f"
		}
#line 236 "slasq5.f"
	    } else {
#line 237 "slasq5.f"
		i__1 = *n0 - 3 << 2;
#line 237 "slasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 238 "slasq5.f"
		    z__[j4 - 3] = d__ + z__[j4];
#line 239 "slasq5.f"
		    if (d__ < 0.) {
#line 240 "slasq5.f"
			return 0;
#line 241 "slasq5.f"
		    } else {
#line 242 "slasq5.f"
			z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
#line 243 "slasq5.f"
			d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
#line 244 "slasq5.f"
		    }
#line 245 "slasq5.f"
		    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 246 "slasq5.f"
		    d__1 = emin, d__2 = z__[j4 - 1];
#line 246 "slasq5.f"
		    emin = min(d__1,d__2);
#line 247 "slasq5.f"
/* L40: */
#line 247 "slasq5.f"
		}
#line 248 "slasq5.f"
	    }

/*     Unroll last two steps. */

#line 252 "slasq5.f"
	    *dnm2 = d__;
#line 253 "slasq5.f"
	    *dmin2 = *dmin__;
#line 254 "slasq5.f"
	    j4 = (*n0 - 2 << 2) - *pp;
#line 255 "slasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 256 "slasq5.f"
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 257 "slasq5.f"
	    if (*dnm2 < 0.) {
#line 258 "slasq5.f"
		return 0;
#line 259 "slasq5.f"
	    } else {
#line 260 "slasq5.f"
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 261 "slasq5.f"
		*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
#line 262 "slasq5.f"
	    }
#line 263 "slasq5.f"
	    *dmin__ = min(*dmin__,*dnm1);

#line 265 "slasq5.f"
	    *dmin1 = *dmin__;
#line 266 "slasq5.f"
	    j4 += 4;
#line 267 "slasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 268 "slasq5.f"
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 269 "slasq5.f"
	    if (*dnm1 < 0.) {
#line 270 "slasq5.f"
		return 0;
#line 271 "slasq5.f"
	    } else {
#line 272 "slasq5.f"
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 273 "slasq5.f"
		*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
#line 274 "slasq5.f"
	    }
#line 275 "slasq5.f"
	    *dmin__ = min(*dmin__,*dn);

#line 277 "slasq5.f"
	}

#line 279 "slasq5.f"
    } else {
/*     This is the version that sets d's to zero if they are small enough */
#line 281 "slasq5.f"
	j4 = (*i0 << 2) + *pp - 3;
#line 282 "slasq5.f"
	emin = z__[j4 + 4];
#line 283 "slasq5.f"
	d__ = z__[j4] - *tau;
#line 284 "slasq5.f"
	*dmin__ = d__;
#line 285 "slasq5.f"
	*dmin1 = -z__[j4];
#line 286 "slasq5.f"
	if (*ieee) {

/*     Code for IEEE arithmetic. */

#line 290 "slasq5.f"
	    if (*pp == 0) {
#line 291 "slasq5.f"
		i__1 = *n0 - 3 << 2;
#line 291 "slasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 292 "slasq5.f"
		    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 293 "slasq5.f"
		    temp = z__[j4 + 1] / z__[j4 - 2];
#line 294 "slasq5.f"
		    d__ = d__ * temp - *tau;
#line 295 "slasq5.f"
		    if (d__ < dthresh) {
#line 295 "slasq5.f"
			d__ = 0.;
#line 295 "slasq5.f"
		    }
#line 296 "slasq5.f"
		    *dmin__ = min(*dmin__,d__);
#line 297 "slasq5.f"
		    z__[j4] = z__[j4 - 1] * temp;
/* Computing MIN */
#line 298 "slasq5.f"
		    d__1 = z__[j4];
#line 298 "slasq5.f"
		    emin = min(d__1,emin);
#line 299 "slasq5.f"
/* L50: */
#line 299 "slasq5.f"
		}
#line 300 "slasq5.f"
	    } else {
#line 301 "slasq5.f"
		i__1 = *n0 - 3 << 2;
#line 301 "slasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 302 "slasq5.f"
		    z__[j4 - 3] = d__ + z__[j4];
#line 303 "slasq5.f"
		    temp = z__[j4 + 2] / z__[j4 - 3];
#line 304 "slasq5.f"
		    d__ = d__ * temp - *tau;
#line 305 "slasq5.f"
		    if (d__ < dthresh) {
#line 305 "slasq5.f"
			d__ = 0.;
#line 305 "slasq5.f"
		    }
#line 306 "slasq5.f"
		    *dmin__ = min(*dmin__,d__);
#line 307 "slasq5.f"
		    z__[j4 - 1] = z__[j4] * temp;
/* Computing MIN */
#line 308 "slasq5.f"
		    d__1 = z__[j4 - 1];
#line 308 "slasq5.f"
		    emin = min(d__1,emin);
#line 309 "slasq5.f"
/* L60: */
#line 309 "slasq5.f"
		}
#line 310 "slasq5.f"
	    }

/*     Unroll last two steps. */

#line 314 "slasq5.f"
	    *dnm2 = d__;
#line 315 "slasq5.f"
	    *dmin2 = *dmin__;
#line 316 "slasq5.f"
	    j4 = (*n0 - 2 << 2) - *pp;
#line 317 "slasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 318 "slasq5.f"
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 319 "slasq5.f"
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 320 "slasq5.f"
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
#line 321 "slasq5.f"
	    *dmin__ = min(*dmin__,*dnm1);

#line 323 "slasq5.f"
	    *dmin1 = *dmin__;
#line 324 "slasq5.f"
	    j4 += 4;
#line 325 "slasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 326 "slasq5.f"
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 327 "slasq5.f"
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 328 "slasq5.f"
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
#line 329 "slasq5.f"
	    *dmin__ = min(*dmin__,*dn);

#line 331 "slasq5.f"
	} else {

/*     Code for non IEEE arithmetic. */

#line 335 "slasq5.f"
	    if (*pp == 0) {
#line 336 "slasq5.f"
		i__1 = *n0 - 3 << 2;
#line 336 "slasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 337 "slasq5.f"
		    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 338 "slasq5.f"
		    if (d__ < 0.) {
#line 339 "slasq5.f"
			return 0;
#line 340 "slasq5.f"
		    } else {
#line 341 "slasq5.f"
			z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
#line 342 "slasq5.f"
			d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
#line 343 "slasq5.f"
		    }
#line 344 "slasq5.f"
		    if (d__ < dthresh) {
#line 344 "slasq5.f"
			d__ = 0.;
#line 344 "slasq5.f"
		    }
#line 345 "slasq5.f"
		    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 346 "slasq5.f"
		    d__1 = emin, d__2 = z__[j4];
#line 346 "slasq5.f"
		    emin = min(d__1,d__2);
#line 347 "slasq5.f"
/* L70: */
#line 347 "slasq5.f"
		}
#line 348 "slasq5.f"
	    } else {
#line 349 "slasq5.f"
		i__1 = *n0 - 3 << 2;
#line 349 "slasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 350 "slasq5.f"
		    z__[j4 - 3] = d__ + z__[j4];
#line 351 "slasq5.f"
		    if (d__ < 0.) {
#line 352 "slasq5.f"
			return 0;
#line 353 "slasq5.f"
		    } else {
#line 354 "slasq5.f"
			z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
#line 355 "slasq5.f"
			d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
#line 356 "slasq5.f"
		    }
#line 357 "slasq5.f"
		    if (d__ < dthresh) {
#line 357 "slasq5.f"
			d__ = 0.;
#line 357 "slasq5.f"
		    }
#line 358 "slasq5.f"
		    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 359 "slasq5.f"
		    d__1 = emin, d__2 = z__[j4 - 1];
#line 359 "slasq5.f"
		    emin = min(d__1,d__2);
#line 360 "slasq5.f"
/* L80: */
#line 360 "slasq5.f"
		}
#line 361 "slasq5.f"
	    }

/*     Unroll last two steps. */

#line 365 "slasq5.f"
	    *dnm2 = d__;
#line 366 "slasq5.f"
	    *dmin2 = *dmin__;
#line 367 "slasq5.f"
	    j4 = (*n0 - 2 << 2) - *pp;
#line 368 "slasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 369 "slasq5.f"
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 370 "slasq5.f"
	    if (*dnm2 < 0.) {
#line 371 "slasq5.f"
		return 0;
#line 372 "slasq5.f"
	    } else {
#line 373 "slasq5.f"
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 374 "slasq5.f"
		*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
#line 375 "slasq5.f"
	    }
#line 376 "slasq5.f"
	    *dmin__ = min(*dmin__,*dnm1);

#line 378 "slasq5.f"
	    *dmin1 = *dmin__;
#line 379 "slasq5.f"
	    j4 += 4;
#line 380 "slasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 381 "slasq5.f"
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 382 "slasq5.f"
	    if (*dnm1 < 0.) {
#line 383 "slasq5.f"
		return 0;
#line 384 "slasq5.f"
	    } else {
#line 385 "slasq5.f"
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 386 "slasq5.f"
		*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
#line 387 "slasq5.f"
	    }
#line 388 "slasq5.f"
	    *dmin__ = min(*dmin__,*dn);

#line 390 "slasq5.f"
	}

#line 392 "slasq5.f"
    }
#line 393 "slasq5.f"
    z__[j4 + 2] = *dn;
#line 394 "slasq5.f"
    z__[(*n0 << 2) - *pp] = emin;
#line 395 "slasq5.f"
    return 0;

/*     End of SLASQ5 */

} /* slasq5_ */

