#line 1 "slaic1.f"
/* slaic1.f -- translated by f2c (version 20100827).
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

#line 1 "slaic1.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = 1.;

/* > \brief \b SLAIC1 applies one step of incremental condition estimation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAIC1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaic1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaic1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaic1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            J, JOB */
/*       REAL               C, GAMMA, S, SEST, SESTPR */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               W( J ), X( J ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAIC1 applies one step of incremental condition estimation in */
/* > its simplest version: */
/* > */
/* > Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/* > lower triangular matrix L, such that */
/* >          twonorm(L*x) = sest */
/* > Then SLAIC1 computes sestpr, s, c such that */
/* > the vector */
/* >                 [ s*x ] */
/* >          xhat = [  c  ] */
/* > is an approximate singular vector of */
/* >                 [ L      0  ] */
/* >          Lhat = [ w**T gamma ] */
/* > in the sense that */
/* >          twonorm(Lhat*xhat) = sestpr. */
/* > */
/* > Depending on JOB, an estimate for the largest or smallest singular */
/* > value is computed. */
/* > */
/* > Note that [s c]**T and sestpr**2 is an eigenpair of the system */
/* > */
/* >     diag(sest*sest, 0) + [alpha  gamma] * [ alpha ] */
/* >                                           [ gamma ] */
/* > */
/* > where  alpha =  x**T*w. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is INTEGER */
/* >          = 1: an estimate for the largest singular value is computed. */
/* >          = 2: an estimate for the smallest singular value is computed. */
/* > \endverbatim */
/* > */
/* > \param[in] J */
/* > \verbatim */
/* >          J is INTEGER */
/* >          Length of X and W */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array, dimension (J) */
/* >          The j-vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] SEST */
/* > \verbatim */
/* >          SEST is REAL */
/* >          Estimated singular value of j by j matrix L */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is REAL array, dimension (J) */
/* >          The j-vector w. */
/* > \endverbatim */
/* > */
/* > \param[in] GAMMA */
/* > \verbatim */
/* >          GAMMA is REAL */
/* >          The diagonal element gamma. */
/* > \endverbatim */
/* > */
/* > \param[out] SESTPR */
/* > \verbatim */
/* >          SESTPR is REAL */
/* >          Estimated singular value of (j+1) by (j+1) matrix Lhat. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL */
/* >          Sine needed in forming xhat. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is REAL */
/* >          Cosine needed in forming xhat. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slaic1_(integer *job, integer *j, doublereal *x, 
	doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
	sestpr, doublereal *s, doublereal *c__)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal b, t, s1, s2, eps, tmp, sine;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal test, zeta1, zeta2, alpha, norma, absgam, absalp;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal cosine, absest;


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
/*     .. Executable Statements .. */

#line 171 "slaic1.f"
    /* Parameter adjustments */
#line 171 "slaic1.f"
    --w;
#line 171 "slaic1.f"
    --x;
#line 171 "slaic1.f"

#line 171 "slaic1.f"
    /* Function Body */
#line 171 "slaic1.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 172 "slaic1.f"
    alpha = sdot_(j, &x[1], &c__1, &w[1], &c__1);

#line 174 "slaic1.f"
    absalp = abs(alpha);
#line 175 "slaic1.f"
    absgam = abs(*gamma);
#line 176 "slaic1.f"
    absest = abs(*sest);

#line 178 "slaic1.f"
    if (*job == 1) {

/*        Estimating largest singular value */

/*        special cases */

#line 184 "slaic1.f"
	if (*sest == 0.) {
#line 185 "slaic1.f"
	    s1 = max(absgam,absalp);
#line 186 "slaic1.f"
	    if (s1 == 0.) {
#line 187 "slaic1.f"
		*s = 0.;
#line 188 "slaic1.f"
		*c__ = 1.;
#line 189 "slaic1.f"
		*sestpr = 0.;
#line 190 "slaic1.f"
	    } else {
#line 191 "slaic1.f"
		*s = alpha / s1;
#line 192 "slaic1.f"
		*c__ = *gamma / s1;
#line 193 "slaic1.f"
		tmp = sqrt(*s * *s + *c__ * *c__);
#line 194 "slaic1.f"
		*s /= tmp;
#line 195 "slaic1.f"
		*c__ /= tmp;
#line 196 "slaic1.f"
		*sestpr = s1 * tmp;
#line 197 "slaic1.f"
	    }
#line 198 "slaic1.f"
	    return 0;
#line 199 "slaic1.f"
	} else if (absgam <= eps * absest) {
#line 200 "slaic1.f"
	    *s = 1.;
#line 201 "slaic1.f"
	    *c__ = 0.;
#line 202 "slaic1.f"
	    tmp = max(absest,absalp);
#line 203 "slaic1.f"
	    s1 = absest / tmp;
#line 204 "slaic1.f"
	    s2 = absalp / tmp;
#line 205 "slaic1.f"
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
#line 206 "slaic1.f"
	    return 0;
#line 207 "slaic1.f"
	} else if (absalp <= eps * absest) {
#line 208 "slaic1.f"
	    s1 = absgam;
#line 209 "slaic1.f"
	    s2 = absest;
#line 210 "slaic1.f"
	    if (s1 <= s2) {
#line 211 "slaic1.f"
		*s = 1.;
#line 212 "slaic1.f"
		*c__ = 0.;
#line 213 "slaic1.f"
		*sestpr = s2;
#line 214 "slaic1.f"
	    } else {
#line 215 "slaic1.f"
		*s = 0.;
#line 216 "slaic1.f"
		*c__ = 1.;
#line 217 "slaic1.f"
		*sestpr = s1;
#line 218 "slaic1.f"
	    }
#line 219 "slaic1.f"
	    return 0;
#line 220 "slaic1.f"
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
#line 221 "slaic1.f"
	    s1 = absgam;
#line 222 "slaic1.f"
	    s2 = absalp;
#line 223 "slaic1.f"
	    if (s1 <= s2) {
#line 224 "slaic1.f"
		tmp = s1 / s2;
#line 225 "slaic1.f"
		*s = sqrt(tmp * tmp + 1.);
#line 226 "slaic1.f"
		*sestpr = s2 * *s;
#line 227 "slaic1.f"
		*c__ = *gamma / s2 / *s;
#line 228 "slaic1.f"
		*s = d_sign(&c_b5, &alpha) / *s;
#line 229 "slaic1.f"
	    } else {
#line 230 "slaic1.f"
		tmp = s2 / s1;
#line 231 "slaic1.f"
		*c__ = sqrt(tmp * tmp + 1.);
#line 232 "slaic1.f"
		*sestpr = s1 * *c__;
#line 233 "slaic1.f"
		*s = alpha / s1 / *c__;
#line 234 "slaic1.f"
		*c__ = d_sign(&c_b5, gamma) / *c__;
#line 235 "slaic1.f"
	    }
#line 236 "slaic1.f"
	    return 0;
#line 237 "slaic1.f"
	} else {

/*           normal case */

#line 241 "slaic1.f"
	    zeta1 = alpha / absest;
#line 242 "slaic1.f"
	    zeta2 = *gamma / absest;

#line 244 "slaic1.f"
	    b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
#line 245 "slaic1.f"
	    *c__ = zeta1 * zeta1;
#line 246 "slaic1.f"
	    if (b > 0.) {
#line 247 "slaic1.f"
		t = *c__ / (b + sqrt(b * b + *c__));
#line 248 "slaic1.f"
	    } else {
#line 249 "slaic1.f"
		t = sqrt(b * b + *c__) - b;
#line 250 "slaic1.f"
	    }

#line 252 "slaic1.f"
	    sine = -zeta1 / t;
#line 253 "slaic1.f"
	    cosine = -zeta2 / (t + 1.);
#line 254 "slaic1.f"
	    tmp = sqrt(sine * sine + cosine * cosine);
#line 255 "slaic1.f"
	    *s = sine / tmp;
#line 256 "slaic1.f"
	    *c__ = cosine / tmp;
#line 257 "slaic1.f"
	    *sestpr = sqrt(t + 1.) * absest;
#line 258 "slaic1.f"
	    return 0;
#line 259 "slaic1.f"
	}

#line 261 "slaic1.f"
    } else if (*job == 2) {

/*        Estimating smallest singular value */

/*        special cases */

#line 267 "slaic1.f"
	if (*sest == 0.) {
#line 268 "slaic1.f"
	    *sestpr = 0.;
#line 269 "slaic1.f"
	    if (max(absgam,absalp) == 0.) {
#line 270 "slaic1.f"
		sine = 1.;
#line 271 "slaic1.f"
		cosine = 0.;
#line 272 "slaic1.f"
	    } else {
#line 273 "slaic1.f"
		sine = -(*gamma);
#line 274 "slaic1.f"
		cosine = alpha;
#line 275 "slaic1.f"
	    }
/* Computing MAX */
#line 276 "slaic1.f"
	    d__1 = abs(sine), d__2 = abs(cosine);
#line 276 "slaic1.f"
	    s1 = max(d__1,d__2);
#line 277 "slaic1.f"
	    *s = sine / s1;
#line 278 "slaic1.f"
	    *c__ = cosine / s1;
#line 279 "slaic1.f"
	    tmp = sqrt(*s * *s + *c__ * *c__);
#line 280 "slaic1.f"
	    *s /= tmp;
#line 281 "slaic1.f"
	    *c__ /= tmp;
#line 282 "slaic1.f"
	    return 0;
#line 283 "slaic1.f"
	} else if (absgam <= eps * absest) {
#line 284 "slaic1.f"
	    *s = 0.;
#line 285 "slaic1.f"
	    *c__ = 1.;
#line 286 "slaic1.f"
	    *sestpr = absgam;
#line 287 "slaic1.f"
	    return 0;
#line 288 "slaic1.f"
	} else if (absalp <= eps * absest) {
#line 289 "slaic1.f"
	    s1 = absgam;
#line 290 "slaic1.f"
	    s2 = absest;
#line 291 "slaic1.f"
	    if (s1 <= s2) {
#line 292 "slaic1.f"
		*s = 0.;
#line 293 "slaic1.f"
		*c__ = 1.;
#line 294 "slaic1.f"
		*sestpr = s1;
#line 295 "slaic1.f"
	    } else {
#line 296 "slaic1.f"
		*s = 1.;
#line 297 "slaic1.f"
		*c__ = 0.;
#line 298 "slaic1.f"
		*sestpr = s2;
#line 299 "slaic1.f"
	    }
#line 300 "slaic1.f"
	    return 0;
#line 301 "slaic1.f"
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
#line 302 "slaic1.f"
	    s1 = absgam;
#line 303 "slaic1.f"
	    s2 = absalp;
#line 304 "slaic1.f"
	    if (s1 <= s2) {
#line 305 "slaic1.f"
		tmp = s1 / s2;
#line 306 "slaic1.f"
		*c__ = sqrt(tmp * tmp + 1.);
#line 307 "slaic1.f"
		*sestpr = absest * (tmp / *c__);
#line 308 "slaic1.f"
		*s = -(*gamma / s2) / *c__;
#line 309 "slaic1.f"
		*c__ = d_sign(&c_b5, &alpha) / *c__;
#line 310 "slaic1.f"
	    } else {
#line 311 "slaic1.f"
		tmp = s2 / s1;
#line 312 "slaic1.f"
		*s = sqrt(tmp * tmp + 1.);
#line 313 "slaic1.f"
		*sestpr = absest / *s;
#line 314 "slaic1.f"
		*c__ = alpha / s1 / *s;
#line 315 "slaic1.f"
		*s = -d_sign(&c_b5, gamma) / *s;
#line 316 "slaic1.f"
	    }
#line 317 "slaic1.f"
	    return 0;
#line 318 "slaic1.f"
	} else {

/*           normal case */

#line 322 "slaic1.f"
	    zeta1 = alpha / absest;
#line 323 "slaic1.f"
	    zeta2 = *gamma / absest;

/* Computing MAX */
#line 325 "slaic1.f"
	    d__3 = zeta1 * zeta1 + 1. + (d__1 = zeta1 * zeta2, abs(d__1)), 
		    d__4 = (d__2 = zeta1 * zeta2, abs(d__2)) + zeta2 * zeta2;
#line 325 "slaic1.f"
	    norma = max(d__3,d__4);

/*           See if root is closer to zero or to ONE */

#line 330 "slaic1.f"
	    test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
#line 331 "slaic1.f"
	    if (test >= 0.) {

/*              root is close to zero, compute directly */

#line 335 "slaic1.f"
		b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
#line 336 "slaic1.f"
		*c__ = zeta2 * zeta2;
#line 337 "slaic1.f"
		t = *c__ / (b + sqrt((d__1 = b * b - *c__, abs(d__1))));
#line 338 "slaic1.f"
		sine = zeta1 / (1. - t);
#line 339 "slaic1.f"
		cosine = -zeta2 / t;
#line 340 "slaic1.f"
		*sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
#line 341 "slaic1.f"
	    } else {

/*              root is closer to ONE, shift by that amount */

#line 345 "slaic1.f"
		b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
#line 346 "slaic1.f"
		*c__ = zeta1 * zeta1;
#line 347 "slaic1.f"
		if (b >= 0.) {
#line 348 "slaic1.f"
		    t = -(*c__) / (b + sqrt(b * b + *c__));
#line 349 "slaic1.f"
		} else {
#line 350 "slaic1.f"
		    t = b - sqrt(b * b + *c__);
#line 351 "slaic1.f"
		}
#line 352 "slaic1.f"
		sine = -zeta1 / t;
#line 353 "slaic1.f"
		cosine = -zeta2 / (t + 1.);
#line 354 "slaic1.f"
		*sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
#line 355 "slaic1.f"
	    }
#line 356 "slaic1.f"
	    tmp = sqrt(sine * sine + cosine * cosine);
#line 357 "slaic1.f"
	    *s = sine / tmp;
#line 358 "slaic1.f"
	    *c__ = cosine / tmp;
#line 359 "slaic1.f"
	    return 0;

#line 361 "slaic1.f"
	}
#line 362 "slaic1.f"
    }
#line 363 "slaic1.f"
    return 0;

/*     End of SLAIC1 */

} /* slaic1_ */

