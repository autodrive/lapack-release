#line 1 "dlasq5.f"
/* dlasq5.f -- translated by f2c (version 20100827).
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

#line 1 "dlasq5.f"
/* > \brief \b DLASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASQ5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, */
/*                          DNM1, DNM2, IEEE, EPS ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            IEEE */
/*       INTEGER            I0, N0, PP */
/*       DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, SIGMA, EPS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ5 computes one dqds transform in ping-pong form, one */
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
/* >          Z is DOUBLE PRECISION array, dimension ( 4*N ) */
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
/* >          TAU is DOUBLE PRECISION */
/* >        This is the shift. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >        This is the accumulated shift up to this step. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN */
/* > \verbatim */
/* >          DMIN is DOUBLE PRECISION */
/* >        Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN1 */
/* > \verbatim */
/* >          DMIN1 is DOUBLE PRECISION */
/* >        Minimum value of d, excluding D( N0 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN2 */
/* > \verbatim */
/* >          DMIN2 is DOUBLE PRECISION */
/* >        Minimum value of d, excluding D( N0 ) and D( N0-1 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DN */
/* > \verbatim */
/* >          DN is DOUBLE PRECISION */
/* >        d(N0), the last value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DNM1 */
/* > \verbatim */
/* >          DNM1 is DOUBLE PRECISION */
/* >        d(N0-1). */
/* > \endverbatim */
/* > */
/* > \param[out] DNM2 */
/* > \verbatim */
/* >          DNM2 is DOUBLE PRECISION */
/* >        d(N0-2). */
/* > \endverbatim */
/* > */
/* > \param[in] IEEE */
/* > \verbatim */
/* >          IEEE is LOGICAL */
/* >        Flag for IEEE or non IEEE arithmetic. */
/* > \endverbatim */

/* > \param[in] EPS */
/* > \verbatim */
/* >          EPS is DOUBLE PRECISION */
/* >        This is the value of epsilon used. */
/* > \endverbatim */
/* > */
/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlasq5_(integer *i0, integer *n0, doublereal *z__, 
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


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 177 "dlasq5.f"
    /* Parameter adjustments */
#line 177 "dlasq5.f"
    --z__;
#line 177 "dlasq5.f"

#line 177 "dlasq5.f"
    /* Function Body */
#line 177 "dlasq5.f"
    if (*n0 - *i0 - 1 <= 0) {
#line 177 "dlasq5.f"
	return 0;
#line 177 "dlasq5.f"
    }

#line 180 "dlasq5.f"
    dthresh = *eps * (*sigma + *tau);
#line 181 "dlasq5.f"
    if (*tau < dthresh * .5) {
#line 181 "dlasq5.f"
	*tau = 0.;
#line 181 "dlasq5.f"
    }
#line 182 "dlasq5.f"
    if (*tau != 0.) {
#line 183 "dlasq5.f"
	j4 = (*i0 << 2) + *pp - 3;
#line 184 "dlasq5.f"
	emin = z__[j4 + 4];
#line 185 "dlasq5.f"
	d__ = z__[j4] - *tau;
#line 186 "dlasq5.f"
	*dmin__ = d__;
#line 187 "dlasq5.f"
	*dmin1 = -z__[j4];

#line 189 "dlasq5.f"
	if (*ieee) {

/*        Code for IEEE arithmetic. */

#line 193 "dlasq5.f"
	    if (*pp == 0) {
#line 194 "dlasq5.f"
		i__1 = *n0 - 3 << 2;
#line 194 "dlasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 195 "dlasq5.f"
		    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 196 "dlasq5.f"
		    temp = z__[j4 + 1] / z__[j4 - 2];
#line 197 "dlasq5.f"
		    d__ = d__ * temp - *tau;
#line 198 "dlasq5.f"
		    *dmin__ = min(*dmin__,d__);
#line 199 "dlasq5.f"
		    z__[j4] = z__[j4 - 1] * temp;
/* Computing MIN */
#line 200 "dlasq5.f"
		    d__1 = z__[j4];
#line 200 "dlasq5.f"
		    emin = min(d__1,emin);
#line 201 "dlasq5.f"
/* L10: */
#line 201 "dlasq5.f"
		}
#line 202 "dlasq5.f"
	    } else {
#line 203 "dlasq5.f"
		i__1 = *n0 - 3 << 2;
#line 203 "dlasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 204 "dlasq5.f"
		    z__[j4 - 3] = d__ + z__[j4];
#line 205 "dlasq5.f"
		    temp = z__[j4 + 2] / z__[j4 - 3];
#line 206 "dlasq5.f"
		    d__ = d__ * temp - *tau;
#line 207 "dlasq5.f"
		    *dmin__ = min(*dmin__,d__);
#line 208 "dlasq5.f"
		    z__[j4 - 1] = z__[j4] * temp;
/* Computing MIN */
#line 209 "dlasq5.f"
		    d__1 = z__[j4 - 1];
#line 209 "dlasq5.f"
		    emin = min(d__1,emin);
#line 210 "dlasq5.f"
/* L20: */
#line 210 "dlasq5.f"
		}
#line 211 "dlasq5.f"
	    }

/*        Unroll last two steps. */

#line 215 "dlasq5.f"
	    *dnm2 = d__;
#line 216 "dlasq5.f"
	    *dmin2 = *dmin__;
#line 217 "dlasq5.f"
	    j4 = (*n0 - 2 << 2) - *pp;
#line 218 "dlasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 219 "dlasq5.f"
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 220 "dlasq5.f"
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 221 "dlasq5.f"
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
#line 222 "dlasq5.f"
	    *dmin__ = min(*dmin__,*dnm1);

#line 224 "dlasq5.f"
	    *dmin1 = *dmin__;
#line 225 "dlasq5.f"
	    j4 += 4;
#line 226 "dlasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 227 "dlasq5.f"
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 228 "dlasq5.f"
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 229 "dlasq5.f"
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
#line 230 "dlasq5.f"
	    *dmin__ = min(*dmin__,*dn);

#line 232 "dlasq5.f"
	} else {

/*        Code for non IEEE arithmetic. */

#line 236 "dlasq5.f"
	    if (*pp == 0) {
#line 237 "dlasq5.f"
		i__1 = *n0 - 3 << 2;
#line 237 "dlasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 238 "dlasq5.f"
		    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 239 "dlasq5.f"
		    if (d__ < 0.) {
#line 240 "dlasq5.f"
			return 0;
#line 241 "dlasq5.f"
		    } else {
#line 242 "dlasq5.f"
			z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
#line 243 "dlasq5.f"
			d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
#line 244 "dlasq5.f"
		    }
#line 245 "dlasq5.f"
		    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 246 "dlasq5.f"
		    d__1 = emin, d__2 = z__[j4];
#line 246 "dlasq5.f"
		    emin = min(d__1,d__2);
#line 247 "dlasq5.f"
/* L30: */
#line 247 "dlasq5.f"
		}
#line 248 "dlasq5.f"
	    } else {
#line 249 "dlasq5.f"
		i__1 = *n0 - 3 << 2;
#line 249 "dlasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 250 "dlasq5.f"
		    z__[j4 - 3] = d__ + z__[j4];
#line 251 "dlasq5.f"
		    if (d__ < 0.) {
#line 252 "dlasq5.f"
			return 0;
#line 253 "dlasq5.f"
		    } else {
#line 254 "dlasq5.f"
			z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
#line 255 "dlasq5.f"
			d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
#line 256 "dlasq5.f"
		    }
#line 257 "dlasq5.f"
		    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 258 "dlasq5.f"
		    d__1 = emin, d__2 = z__[j4 - 1];
#line 258 "dlasq5.f"
		    emin = min(d__1,d__2);
#line 259 "dlasq5.f"
/* L40: */
#line 259 "dlasq5.f"
		}
#line 260 "dlasq5.f"
	    }

/*        Unroll last two steps. */

#line 264 "dlasq5.f"
	    *dnm2 = d__;
#line 265 "dlasq5.f"
	    *dmin2 = *dmin__;
#line 266 "dlasq5.f"
	    j4 = (*n0 - 2 << 2) - *pp;
#line 267 "dlasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 268 "dlasq5.f"
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 269 "dlasq5.f"
	    if (*dnm2 < 0.) {
#line 270 "dlasq5.f"
		return 0;
#line 271 "dlasq5.f"
	    } else {
#line 272 "dlasq5.f"
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 273 "dlasq5.f"
		*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
#line 274 "dlasq5.f"
	    }
#line 275 "dlasq5.f"
	    *dmin__ = min(*dmin__,*dnm1);

#line 277 "dlasq5.f"
	    *dmin1 = *dmin__;
#line 278 "dlasq5.f"
	    j4 += 4;
#line 279 "dlasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 280 "dlasq5.f"
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 281 "dlasq5.f"
	    if (*dnm1 < 0.) {
#line 282 "dlasq5.f"
		return 0;
#line 283 "dlasq5.f"
	    } else {
#line 284 "dlasq5.f"
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 285 "dlasq5.f"
		*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
#line 286 "dlasq5.f"
	    }
#line 287 "dlasq5.f"
	    *dmin__ = min(*dmin__,*dn);

#line 289 "dlasq5.f"
	}
#line 290 "dlasq5.f"
    } else {
/*     This is the version that sets d's to zero if they are small enough */
#line 292 "dlasq5.f"
	j4 = (*i0 << 2) + *pp - 3;
#line 293 "dlasq5.f"
	emin = z__[j4 + 4];
#line 294 "dlasq5.f"
	d__ = z__[j4] - *tau;
#line 295 "dlasq5.f"
	*dmin__ = d__;
#line 296 "dlasq5.f"
	*dmin1 = -z__[j4];
#line 297 "dlasq5.f"
	if (*ieee) {

/*     Code for IEEE arithmetic. */

#line 301 "dlasq5.f"
	    if (*pp == 0) {
#line 302 "dlasq5.f"
		i__1 = *n0 - 3 << 2;
#line 302 "dlasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 303 "dlasq5.f"
		    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 304 "dlasq5.f"
		    temp = z__[j4 + 1] / z__[j4 - 2];
#line 305 "dlasq5.f"
		    d__ = d__ * temp - *tau;
#line 306 "dlasq5.f"
		    if (d__ < dthresh) {
#line 306 "dlasq5.f"
			d__ = 0.;
#line 306 "dlasq5.f"
		    }
#line 307 "dlasq5.f"
		    *dmin__ = min(*dmin__,d__);
#line 308 "dlasq5.f"
		    z__[j4] = z__[j4 - 1] * temp;
/* Computing MIN */
#line 309 "dlasq5.f"
		    d__1 = z__[j4];
#line 309 "dlasq5.f"
		    emin = min(d__1,emin);
#line 310 "dlasq5.f"
/* L50: */
#line 310 "dlasq5.f"
		}
#line 311 "dlasq5.f"
	    } else {
#line 312 "dlasq5.f"
		i__1 = *n0 - 3 << 2;
#line 312 "dlasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 313 "dlasq5.f"
		    z__[j4 - 3] = d__ + z__[j4];
#line 314 "dlasq5.f"
		    temp = z__[j4 + 2] / z__[j4 - 3];
#line 315 "dlasq5.f"
		    d__ = d__ * temp - *tau;
#line 316 "dlasq5.f"
		    if (d__ < dthresh) {
#line 316 "dlasq5.f"
			d__ = 0.;
#line 316 "dlasq5.f"
		    }
#line 317 "dlasq5.f"
		    *dmin__ = min(*dmin__,d__);
#line 318 "dlasq5.f"
		    z__[j4 - 1] = z__[j4] * temp;
/* Computing MIN */
#line 319 "dlasq5.f"
		    d__1 = z__[j4 - 1];
#line 319 "dlasq5.f"
		    emin = min(d__1,emin);
#line 320 "dlasq5.f"
/* L60: */
#line 320 "dlasq5.f"
		}
#line 321 "dlasq5.f"
	    }

/*     Unroll last two steps. */

#line 325 "dlasq5.f"
	    *dnm2 = d__;
#line 326 "dlasq5.f"
	    *dmin2 = *dmin__;
#line 327 "dlasq5.f"
	    j4 = (*n0 - 2 << 2) - *pp;
#line 328 "dlasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 329 "dlasq5.f"
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 330 "dlasq5.f"
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 331 "dlasq5.f"
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
#line 332 "dlasq5.f"
	    *dmin__ = min(*dmin__,*dnm1);

#line 334 "dlasq5.f"
	    *dmin1 = *dmin__;
#line 335 "dlasq5.f"
	    j4 += 4;
#line 336 "dlasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 337 "dlasq5.f"
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 338 "dlasq5.f"
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 339 "dlasq5.f"
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
#line 340 "dlasq5.f"
	    *dmin__ = min(*dmin__,*dn);

#line 342 "dlasq5.f"
	} else {

/*     Code for non IEEE arithmetic. */

#line 346 "dlasq5.f"
	    if (*pp == 0) {
#line 347 "dlasq5.f"
		i__1 = *n0 - 3 << 2;
#line 347 "dlasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 348 "dlasq5.f"
		    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 349 "dlasq5.f"
		    if (d__ < 0.) {
#line 350 "dlasq5.f"
			return 0;
#line 351 "dlasq5.f"
		    } else {
#line 352 "dlasq5.f"
			z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
#line 353 "dlasq5.f"
			d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
#line 354 "dlasq5.f"
		    }
#line 355 "dlasq5.f"
		    if (d__ < dthresh) {
#line 355 "dlasq5.f"
			d__ = 0.;
#line 355 "dlasq5.f"
		    }
#line 356 "dlasq5.f"
		    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 357 "dlasq5.f"
		    d__1 = emin, d__2 = z__[j4];
#line 357 "dlasq5.f"
		    emin = min(d__1,d__2);
#line 358 "dlasq5.f"
/* L70: */
#line 358 "dlasq5.f"
		}
#line 359 "dlasq5.f"
	    } else {
#line 360 "dlasq5.f"
		i__1 = *n0 - 3 << 2;
#line 360 "dlasq5.f"
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 361 "dlasq5.f"
		    z__[j4 - 3] = d__ + z__[j4];
#line 362 "dlasq5.f"
		    if (d__ < 0.) {
#line 363 "dlasq5.f"
			return 0;
#line 364 "dlasq5.f"
		    } else {
#line 365 "dlasq5.f"
			z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
#line 366 "dlasq5.f"
			d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
#line 367 "dlasq5.f"
		    }
#line 368 "dlasq5.f"
		    if (d__ < dthresh) {
#line 368 "dlasq5.f"
			d__ = 0.;
#line 368 "dlasq5.f"
		    }
#line 369 "dlasq5.f"
		    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 370 "dlasq5.f"
		    d__1 = emin, d__2 = z__[j4 - 1];
#line 370 "dlasq5.f"
		    emin = min(d__1,d__2);
#line 371 "dlasq5.f"
/* L80: */
#line 371 "dlasq5.f"
		}
#line 372 "dlasq5.f"
	    }

/*     Unroll last two steps. */

#line 376 "dlasq5.f"
	    *dnm2 = d__;
#line 377 "dlasq5.f"
	    *dmin2 = *dmin__;
#line 378 "dlasq5.f"
	    j4 = (*n0 - 2 << 2) - *pp;
#line 379 "dlasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 380 "dlasq5.f"
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 381 "dlasq5.f"
	    if (*dnm2 < 0.) {
#line 382 "dlasq5.f"
		return 0;
#line 383 "dlasq5.f"
	    } else {
#line 384 "dlasq5.f"
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 385 "dlasq5.f"
		*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
#line 386 "dlasq5.f"
	    }
#line 387 "dlasq5.f"
	    *dmin__ = min(*dmin__,*dnm1);

#line 389 "dlasq5.f"
	    *dmin1 = *dmin__;
#line 390 "dlasq5.f"
	    j4 += 4;
#line 391 "dlasq5.f"
	    j4p2 = j4 + (*pp << 1) - 1;
#line 392 "dlasq5.f"
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 393 "dlasq5.f"
	    if (*dnm1 < 0.) {
#line 394 "dlasq5.f"
		return 0;
#line 395 "dlasq5.f"
	    } else {
#line 396 "dlasq5.f"
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 397 "dlasq5.f"
		*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
#line 398 "dlasq5.f"
	    }
#line 399 "dlasq5.f"
	    *dmin__ = min(*dmin__,*dn);

#line 401 "dlasq5.f"
	}
#line 402 "dlasq5.f"
    }

#line 404 "dlasq5.f"
    z__[j4 + 2] = *dn;
#line 405 "dlasq5.f"
    z__[(*n0 << 2) - *pp] = emin;
#line 406 "dlasq5.f"
    return 0;

/*     End of DLASQ5 */

} /* dlasq5_ */

