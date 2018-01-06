#line 1 "dlaic1.f"
/* dlaic1.f -- translated by f2c (version 20100827).
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

#line 1 "dlaic1.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = 1.;

/* > \brief \b DLAIC1 applies one step of incremental condition estimation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAIC1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaic1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaic1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaic1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            J, JOB */
/*       DOUBLE PRECISION   C, GAMMA, S, SEST, SESTPR */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   W( J ), X( J ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAIC1 applies one step of incremental condition estimation in */
/* > its simplest version: */
/* > */
/* > Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/* > lower triangular matrix L, such that */
/* >          twonorm(L*x) = sest */
/* > Then DLAIC1 computes sestpr, s, c such that */
/* > the vector */
/* >                 [ s*x ] */
/* >          xhat = [  c  ] */
/* > is an approximate singular vector of */
/* >                 [ L       0  ] */
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
/* >          X is DOUBLE PRECISION array, dimension (J) */
/* >          The j-vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] SEST */
/* > \verbatim */
/* >          SEST is DOUBLE PRECISION */
/* >          Estimated singular value of j by j matrix L */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (J) */
/* >          The j-vector w. */
/* > \endverbatim */
/* > */
/* > \param[in] GAMMA */
/* > \verbatim */
/* >          GAMMA is DOUBLE PRECISION */
/* >          The diagonal element gamma. */
/* > \endverbatim */
/* > */
/* > \param[out] SESTPR */
/* > \verbatim */
/* >          SESTPR is DOUBLE PRECISION */
/* >          Estimated singular value of (j+1) by (j+1) matrix Lhat. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION */
/* >          Sine needed in forming xhat. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* >          Cosine needed in forming xhat. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaic1_(integer *job, integer *j, doublereal *x, 
	doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
	sestpr, doublereal *s, doublereal *c__)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal b, t, s1, s2, eps, tmp;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal sine, test, zeta1, zeta2, alpha, norma;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal absgam, absalp, cosine, absest;


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

#line 171 "dlaic1.f"
    /* Parameter adjustments */
#line 171 "dlaic1.f"
    --w;
#line 171 "dlaic1.f"
    --x;
#line 171 "dlaic1.f"

#line 171 "dlaic1.f"
    /* Function Body */
#line 171 "dlaic1.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 172 "dlaic1.f"
    alpha = ddot_(j, &x[1], &c__1, &w[1], &c__1);

#line 174 "dlaic1.f"
    absalp = abs(alpha);
#line 175 "dlaic1.f"
    absgam = abs(*gamma);
#line 176 "dlaic1.f"
    absest = abs(*sest);

#line 178 "dlaic1.f"
    if (*job == 1) {

/*        Estimating largest singular value */

/*        special cases */

#line 184 "dlaic1.f"
	if (*sest == 0.) {
#line 185 "dlaic1.f"
	    s1 = max(absgam,absalp);
#line 186 "dlaic1.f"
	    if (s1 == 0.) {
#line 187 "dlaic1.f"
		*s = 0.;
#line 188 "dlaic1.f"
		*c__ = 1.;
#line 189 "dlaic1.f"
		*sestpr = 0.;
#line 190 "dlaic1.f"
	    } else {
#line 191 "dlaic1.f"
		*s = alpha / s1;
#line 192 "dlaic1.f"
		*c__ = *gamma / s1;
#line 193 "dlaic1.f"
		tmp = sqrt(*s * *s + *c__ * *c__);
#line 194 "dlaic1.f"
		*s /= tmp;
#line 195 "dlaic1.f"
		*c__ /= tmp;
#line 196 "dlaic1.f"
		*sestpr = s1 * tmp;
#line 197 "dlaic1.f"
	    }
#line 198 "dlaic1.f"
	    return 0;
#line 199 "dlaic1.f"
	} else if (absgam <= eps * absest) {
#line 200 "dlaic1.f"
	    *s = 1.;
#line 201 "dlaic1.f"
	    *c__ = 0.;
#line 202 "dlaic1.f"
	    tmp = max(absest,absalp);
#line 203 "dlaic1.f"
	    s1 = absest / tmp;
#line 204 "dlaic1.f"
	    s2 = absalp / tmp;
#line 205 "dlaic1.f"
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
#line 206 "dlaic1.f"
	    return 0;
#line 207 "dlaic1.f"
	} else if (absalp <= eps * absest) {
#line 208 "dlaic1.f"
	    s1 = absgam;
#line 209 "dlaic1.f"
	    s2 = absest;
#line 210 "dlaic1.f"
	    if (s1 <= s2) {
#line 211 "dlaic1.f"
		*s = 1.;
#line 212 "dlaic1.f"
		*c__ = 0.;
#line 213 "dlaic1.f"
		*sestpr = s2;
#line 214 "dlaic1.f"
	    } else {
#line 215 "dlaic1.f"
		*s = 0.;
#line 216 "dlaic1.f"
		*c__ = 1.;
#line 217 "dlaic1.f"
		*sestpr = s1;
#line 218 "dlaic1.f"
	    }
#line 219 "dlaic1.f"
	    return 0;
#line 220 "dlaic1.f"
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
#line 221 "dlaic1.f"
	    s1 = absgam;
#line 222 "dlaic1.f"
	    s2 = absalp;
#line 223 "dlaic1.f"
	    if (s1 <= s2) {
#line 224 "dlaic1.f"
		tmp = s1 / s2;
#line 225 "dlaic1.f"
		*s = sqrt(tmp * tmp + 1.);
#line 226 "dlaic1.f"
		*sestpr = s2 * *s;
#line 227 "dlaic1.f"
		*c__ = *gamma / s2 / *s;
#line 228 "dlaic1.f"
		*s = d_sign(&c_b5, &alpha) / *s;
#line 229 "dlaic1.f"
	    } else {
#line 230 "dlaic1.f"
		tmp = s2 / s1;
#line 231 "dlaic1.f"
		*c__ = sqrt(tmp * tmp + 1.);
#line 232 "dlaic1.f"
		*sestpr = s1 * *c__;
#line 233 "dlaic1.f"
		*s = alpha / s1 / *c__;
#line 234 "dlaic1.f"
		*c__ = d_sign(&c_b5, gamma) / *c__;
#line 235 "dlaic1.f"
	    }
#line 236 "dlaic1.f"
	    return 0;
#line 237 "dlaic1.f"
	} else {

/*           normal case */

#line 241 "dlaic1.f"
	    zeta1 = alpha / absest;
#line 242 "dlaic1.f"
	    zeta2 = *gamma / absest;

#line 244 "dlaic1.f"
	    b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
#line 245 "dlaic1.f"
	    *c__ = zeta1 * zeta1;
#line 246 "dlaic1.f"
	    if (b > 0.) {
#line 247 "dlaic1.f"
		t = *c__ / (b + sqrt(b * b + *c__));
#line 248 "dlaic1.f"
	    } else {
#line 249 "dlaic1.f"
		t = sqrt(b * b + *c__) - b;
#line 250 "dlaic1.f"
	    }

#line 252 "dlaic1.f"
	    sine = -zeta1 / t;
#line 253 "dlaic1.f"
	    cosine = -zeta2 / (t + 1.);
#line 254 "dlaic1.f"
	    tmp = sqrt(sine * sine + cosine * cosine);
#line 255 "dlaic1.f"
	    *s = sine / tmp;
#line 256 "dlaic1.f"
	    *c__ = cosine / tmp;
#line 257 "dlaic1.f"
	    *sestpr = sqrt(t + 1.) * absest;
#line 258 "dlaic1.f"
	    return 0;
#line 259 "dlaic1.f"
	}

#line 261 "dlaic1.f"
    } else if (*job == 2) {

/*        Estimating smallest singular value */

/*        special cases */

#line 267 "dlaic1.f"
	if (*sest == 0.) {
#line 268 "dlaic1.f"
	    *sestpr = 0.;
#line 269 "dlaic1.f"
	    if (max(absgam,absalp) == 0.) {
#line 270 "dlaic1.f"
		sine = 1.;
#line 271 "dlaic1.f"
		cosine = 0.;
#line 272 "dlaic1.f"
	    } else {
#line 273 "dlaic1.f"
		sine = -(*gamma);
#line 274 "dlaic1.f"
		cosine = alpha;
#line 275 "dlaic1.f"
	    }
/* Computing MAX */
#line 276 "dlaic1.f"
	    d__1 = abs(sine), d__2 = abs(cosine);
#line 276 "dlaic1.f"
	    s1 = max(d__1,d__2);
#line 277 "dlaic1.f"
	    *s = sine / s1;
#line 278 "dlaic1.f"
	    *c__ = cosine / s1;
#line 279 "dlaic1.f"
	    tmp = sqrt(*s * *s + *c__ * *c__);
#line 280 "dlaic1.f"
	    *s /= tmp;
#line 281 "dlaic1.f"
	    *c__ /= tmp;
#line 282 "dlaic1.f"
	    return 0;
#line 283 "dlaic1.f"
	} else if (absgam <= eps * absest) {
#line 284 "dlaic1.f"
	    *s = 0.;
#line 285 "dlaic1.f"
	    *c__ = 1.;
#line 286 "dlaic1.f"
	    *sestpr = absgam;
#line 287 "dlaic1.f"
	    return 0;
#line 288 "dlaic1.f"
	} else if (absalp <= eps * absest) {
#line 289 "dlaic1.f"
	    s1 = absgam;
#line 290 "dlaic1.f"
	    s2 = absest;
#line 291 "dlaic1.f"
	    if (s1 <= s2) {
#line 292 "dlaic1.f"
		*s = 0.;
#line 293 "dlaic1.f"
		*c__ = 1.;
#line 294 "dlaic1.f"
		*sestpr = s1;
#line 295 "dlaic1.f"
	    } else {
#line 296 "dlaic1.f"
		*s = 1.;
#line 297 "dlaic1.f"
		*c__ = 0.;
#line 298 "dlaic1.f"
		*sestpr = s2;
#line 299 "dlaic1.f"
	    }
#line 300 "dlaic1.f"
	    return 0;
#line 301 "dlaic1.f"
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
#line 302 "dlaic1.f"
	    s1 = absgam;
#line 303 "dlaic1.f"
	    s2 = absalp;
#line 304 "dlaic1.f"
	    if (s1 <= s2) {
#line 305 "dlaic1.f"
		tmp = s1 / s2;
#line 306 "dlaic1.f"
		*c__ = sqrt(tmp * tmp + 1.);
#line 307 "dlaic1.f"
		*sestpr = absest * (tmp / *c__);
#line 308 "dlaic1.f"
		*s = -(*gamma / s2) / *c__;
#line 309 "dlaic1.f"
		*c__ = d_sign(&c_b5, &alpha) / *c__;
#line 310 "dlaic1.f"
	    } else {
#line 311 "dlaic1.f"
		tmp = s2 / s1;
#line 312 "dlaic1.f"
		*s = sqrt(tmp * tmp + 1.);
#line 313 "dlaic1.f"
		*sestpr = absest / *s;
#line 314 "dlaic1.f"
		*c__ = alpha / s1 / *s;
#line 315 "dlaic1.f"
		*s = -d_sign(&c_b5, gamma) / *s;
#line 316 "dlaic1.f"
	    }
#line 317 "dlaic1.f"
	    return 0;
#line 318 "dlaic1.f"
	} else {

/*           normal case */

#line 322 "dlaic1.f"
	    zeta1 = alpha / absest;
#line 323 "dlaic1.f"
	    zeta2 = *gamma / absest;

/* Computing MAX */
#line 325 "dlaic1.f"
	    d__3 = zeta1 * zeta1 + 1. + (d__1 = zeta1 * zeta2, abs(d__1)), 
		    d__4 = (d__2 = zeta1 * zeta2, abs(d__2)) + zeta2 * zeta2;
#line 325 "dlaic1.f"
	    norma = max(d__3,d__4);

/*           See if root is closer to zero or to ONE */

#line 330 "dlaic1.f"
	    test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
#line 331 "dlaic1.f"
	    if (test >= 0.) {

/*              root is close to zero, compute directly */

#line 335 "dlaic1.f"
		b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
#line 336 "dlaic1.f"
		*c__ = zeta2 * zeta2;
#line 337 "dlaic1.f"
		t = *c__ / (b + sqrt((d__1 = b * b - *c__, abs(d__1))));
#line 338 "dlaic1.f"
		sine = zeta1 / (1. - t);
#line 339 "dlaic1.f"
		cosine = -zeta2 / t;
#line 340 "dlaic1.f"
		*sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
#line 341 "dlaic1.f"
	    } else {

/*              root is closer to ONE, shift by that amount */

#line 345 "dlaic1.f"
		b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
#line 346 "dlaic1.f"
		*c__ = zeta1 * zeta1;
#line 347 "dlaic1.f"
		if (b >= 0.) {
#line 348 "dlaic1.f"
		    t = -(*c__) / (b + sqrt(b * b + *c__));
#line 349 "dlaic1.f"
		} else {
#line 350 "dlaic1.f"
		    t = b - sqrt(b * b + *c__);
#line 351 "dlaic1.f"
		}
#line 352 "dlaic1.f"
		sine = -zeta1 / t;
#line 353 "dlaic1.f"
		cosine = -zeta2 / (t + 1.);
#line 354 "dlaic1.f"
		*sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
#line 355 "dlaic1.f"
	    }
#line 356 "dlaic1.f"
	    tmp = sqrt(sine * sine + cosine * cosine);
#line 357 "dlaic1.f"
	    *s = sine / tmp;
#line 358 "dlaic1.f"
	    *c__ = cosine / tmp;
#line 359 "dlaic1.f"
	    return 0;

#line 361 "dlaic1.f"
	}
#line 362 "dlaic1.f"
    }
#line 363 "dlaic1.f"
    return 0;

/*     End of DLAIC1 */

} /* dlaic1_ */

