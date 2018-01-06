#line 1 "claic1.f"
/* claic1.f -- translated by f2c (version 20100827).
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

#line 1 "claic1.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLAIC1 applies one step of incremental condition estimation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAIC1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claic1.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claic1.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claic1.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            J, JOB */
/*       REAL               SEST, SESTPR */
/*       COMPLEX            C, GAMMA, S */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            W( J ), X( J ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAIC1 applies one step of incremental condition estimation in */
/* > its simplest version: */
/* > */
/* > Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/* > lower triangular matrix L, such that */
/* >          twonorm(L*x) = sest */
/* > Then CLAIC1 computes sestpr, s, c such that */
/* > the vector */
/* >                 [ s*x ] */
/* >          xhat = [  c  ] */
/* > is an approximate singular vector of */
/* >                 [ L      0  ] */
/* >          Lhat = [ w**H gamma ] */
/* > in the sense that */
/* >          twonorm(Lhat*xhat) = sestpr. */
/* > */
/* > Depending on JOB, an estimate for the largest or smallest singular */
/* > value is computed. */
/* > */
/* > Note that [s c]**H and sestpr**2 is an eigenpair of the system */
/* > */
/* >     diag(sest*sest, 0) + [alpha  gamma] * [ conjg(alpha) ] */
/* >                                           [ conjg(gamma) ] */
/* > */
/* > where  alpha =  x**H*w. */
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
/* >          X is COMPLEX array, dimension (J) */
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
/* >          W is COMPLEX array, dimension (J) */
/* >          The j-vector w. */
/* > \endverbatim */
/* > */
/* > \param[in] GAMMA */
/* > \verbatim */
/* >          GAMMA is COMPLEX */
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
/* >          S is COMPLEX */
/* >          Sine needed in forming xhat. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is COMPLEX */
/* >          Cosine needed in forming xhat. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int claic1_(integer *job, integer *j, doublecomplex *x, 
	doublereal *sest, doublecomplex *w, doublecomplex *gamma, doublereal *
	sestpr, doublecomplex *s, doublecomplex *c__)
{
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_sqrt(doublecomplex *, 
	    doublecomplex *);
    double sqrt(doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal b, t, s1, s2, scl, eps, tmp;
    static doublecomplex sine;
    static doublereal test, zeta1, zeta2;
    static doublecomplex alpha;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal norma, absgam, absalp;
    extern doublereal slamch_(char *, ftnlen);
    static doublecomplex cosine;
    static doublereal absest;


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
/*     .. Executable Statements .. */

#line 175 "claic1.f"
    /* Parameter adjustments */
#line 175 "claic1.f"
    --w;
#line 175 "claic1.f"
    --x;
#line 175 "claic1.f"

#line 175 "claic1.f"
    /* Function Body */
#line 175 "claic1.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 176 "claic1.f"
    cdotc_(&z__1, j, &x[1], &c__1, &w[1], &c__1);
#line 176 "claic1.f"
    alpha.r = z__1.r, alpha.i = z__1.i;

#line 178 "claic1.f"
    absalp = z_abs(&alpha);
#line 179 "claic1.f"
    absgam = z_abs(gamma);
#line 180 "claic1.f"
    absest = abs(*sest);

#line 182 "claic1.f"
    if (*job == 1) {

/*        Estimating largest singular value */

/*        special cases */

#line 188 "claic1.f"
	if (*sest == 0.) {
#line 189 "claic1.f"
	    s1 = max(absgam,absalp);
#line 190 "claic1.f"
	    if (s1 == 0.) {
#line 191 "claic1.f"
		s->r = 0., s->i = 0.;
#line 192 "claic1.f"
		c__->r = 1., c__->i = 0.;
#line 193 "claic1.f"
		*sestpr = 0.;
#line 194 "claic1.f"
	    } else {
#line 195 "claic1.f"
		z__1.r = alpha.r / s1, z__1.i = alpha.i / s1;
#line 195 "claic1.f"
		s->r = z__1.r, s->i = z__1.i;
#line 196 "claic1.f"
		z__1.r = gamma->r / s1, z__1.i = gamma->i / s1;
#line 196 "claic1.f"
		c__->r = z__1.r, c__->i = z__1.i;
#line 197 "claic1.f"
		d_cnjg(&z__4, s);
#line 197 "claic1.f"
		z__3.r = s->r * z__4.r - s->i * z__4.i, z__3.i = s->r * 
			z__4.i + s->i * z__4.r;
#line 197 "claic1.f"
		d_cnjg(&z__6, c__);
#line 197 "claic1.f"
		z__5.r = c__->r * z__6.r - c__->i * z__6.i, z__5.i = c__->r * 
			z__6.i + c__->i * z__6.r;
#line 197 "claic1.f"
		z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
#line 197 "claic1.f"
		z_sqrt(&z__1, &z__2);
#line 197 "claic1.f"
		tmp = z__1.r;
#line 198 "claic1.f"
		z__1.r = s->r / tmp, z__1.i = s->i / tmp;
#line 198 "claic1.f"
		s->r = z__1.r, s->i = z__1.i;
#line 199 "claic1.f"
		z__1.r = c__->r / tmp, z__1.i = c__->i / tmp;
#line 199 "claic1.f"
		c__->r = z__1.r, c__->i = z__1.i;
#line 200 "claic1.f"
		*sestpr = s1 * tmp;
#line 201 "claic1.f"
	    }
#line 202 "claic1.f"
	    return 0;
#line 203 "claic1.f"
	} else if (absgam <= eps * absest) {
#line 204 "claic1.f"
	    s->r = 1., s->i = 0.;
#line 205 "claic1.f"
	    c__->r = 0., c__->i = 0.;
#line 206 "claic1.f"
	    tmp = max(absest,absalp);
#line 207 "claic1.f"
	    s1 = absest / tmp;
#line 208 "claic1.f"
	    s2 = absalp / tmp;
#line 209 "claic1.f"
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
#line 210 "claic1.f"
	    return 0;
#line 211 "claic1.f"
	} else if (absalp <= eps * absest) {
#line 212 "claic1.f"
	    s1 = absgam;
#line 213 "claic1.f"
	    s2 = absest;
#line 214 "claic1.f"
	    if (s1 <= s2) {
#line 215 "claic1.f"
		s->r = 1., s->i = 0.;
#line 216 "claic1.f"
		c__->r = 0., c__->i = 0.;
#line 217 "claic1.f"
		*sestpr = s2;
#line 218 "claic1.f"
	    } else {
#line 219 "claic1.f"
		s->r = 0., s->i = 0.;
#line 220 "claic1.f"
		c__->r = 1., c__->i = 0.;
#line 221 "claic1.f"
		*sestpr = s1;
#line 222 "claic1.f"
	    }
#line 223 "claic1.f"
	    return 0;
#line 224 "claic1.f"
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
#line 225 "claic1.f"
	    s1 = absgam;
#line 226 "claic1.f"
	    s2 = absalp;
#line 227 "claic1.f"
	    if (s1 <= s2) {
#line 228 "claic1.f"
		tmp = s1 / s2;
#line 229 "claic1.f"
		scl = sqrt(tmp * tmp + 1.);
#line 230 "claic1.f"
		*sestpr = s2 * scl;
#line 231 "claic1.f"
		z__2.r = alpha.r / s2, z__2.i = alpha.i / s2;
#line 231 "claic1.f"
		z__1.r = z__2.r / scl, z__1.i = z__2.i / scl;
#line 231 "claic1.f"
		s->r = z__1.r, s->i = z__1.i;
#line 232 "claic1.f"
		z__2.r = gamma->r / s2, z__2.i = gamma->i / s2;
#line 232 "claic1.f"
		z__1.r = z__2.r / scl, z__1.i = z__2.i / scl;
#line 232 "claic1.f"
		c__->r = z__1.r, c__->i = z__1.i;
#line 233 "claic1.f"
	    } else {
#line 234 "claic1.f"
		tmp = s2 / s1;
#line 235 "claic1.f"
		scl = sqrt(tmp * tmp + 1.);
#line 236 "claic1.f"
		*sestpr = s1 * scl;
#line 237 "claic1.f"
		z__2.r = alpha.r / s1, z__2.i = alpha.i / s1;
#line 237 "claic1.f"
		z__1.r = z__2.r / scl, z__1.i = z__2.i / scl;
#line 237 "claic1.f"
		s->r = z__1.r, s->i = z__1.i;
#line 238 "claic1.f"
		z__2.r = gamma->r / s1, z__2.i = gamma->i / s1;
#line 238 "claic1.f"
		z__1.r = z__2.r / scl, z__1.i = z__2.i / scl;
#line 238 "claic1.f"
		c__->r = z__1.r, c__->i = z__1.i;
#line 239 "claic1.f"
	    }
#line 240 "claic1.f"
	    return 0;
#line 241 "claic1.f"
	} else {

/*           normal case */

#line 245 "claic1.f"
	    zeta1 = absalp / absest;
#line 246 "claic1.f"
	    zeta2 = absgam / absest;

#line 248 "claic1.f"
	    b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
#line 249 "claic1.f"
	    d__1 = zeta1 * zeta1;
#line 249 "claic1.f"
	    c__->r = d__1, c__->i = 0.;
#line 250 "claic1.f"
	    if (b > 0.) {
#line 251 "claic1.f"
		d__1 = b * b;
#line 251 "claic1.f"
		z__4.r = d__1 + c__->r, z__4.i = c__->i;
#line 251 "claic1.f"
		z_sqrt(&z__3, &z__4);
#line 251 "claic1.f"
		z__2.r = b + z__3.r, z__2.i = z__3.i;
#line 251 "claic1.f"
		z_div(&z__1, c__, &z__2);
#line 251 "claic1.f"
		t = z__1.r;
#line 252 "claic1.f"
	    } else {
#line 253 "claic1.f"
		d__1 = b * b;
#line 253 "claic1.f"
		z__3.r = d__1 + c__->r, z__3.i = c__->i;
#line 253 "claic1.f"
		z_sqrt(&z__2, &z__3);
#line 253 "claic1.f"
		z__1.r = z__2.r - b, z__1.i = z__2.i;
#line 253 "claic1.f"
		t = z__1.r;
#line 254 "claic1.f"
	    }

#line 256 "claic1.f"
	    z__3.r = alpha.r / absest, z__3.i = alpha.i / absest;
#line 256 "claic1.f"
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 256 "claic1.f"
	    z__1.r = z__2.r / t, z__1.i = z__2.i / t;
#line 256 "claic1.f"
	    sine.r = z__1.r, sine.i = z__1.i;
#line 257 "claic1.f"
	    z__3.r = gamma->r / absest, z__3.i = gamma->i / absest;
#line 257 "claic1.f"
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 257 "claic1.f"
	    d__1 = t + 1.;
#line 257 "claic1.f"
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 257 "claic1.f"
	    cosine.r = z__1.r, cosine.i = z__1.i;
#line 258 "claic1.f"
	    d_cnjg(&z__4, &sine);
#line 258 "claic1.f"
	    z__3.r = sine.r * z__4.r - sine.i * z__4.i, z__3.i = sine.r * 
		    z__4.i + sine.i * z__4.r;
#line 258 "claic1.f"
	    d_cnjg(&z__6, &cosine);
#line 258 "claic1.f"
	    z__5.r = cosine.r * z__6.r - cosine.i * z__6.i, z__5.i = cosine.r 
		    * z__6.i + cosine.i * z__6.r;
#line 258 "claic1.f"
	    z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
#line 258 "claic1.f"
	    z_sqrt(&z__1, &z__2);
#line 258 "claic1.f"
	    tmp = z__1.r;
#line 259 "claic1.f"
	    z__1.r = sine.r / tmp, z__1.i = sine.i / tmp;
#line 259 "claic1.f"
	    s->r = z__1.r, s->i = z__1.i;
#line 260 "claic1.f"
	    z__1.r = cosine.r / tmp, z__1.i = cosine.i / tmp;
#line 260 "claic1.f"
	    c__->r = z__1.r, c__->i = z__1.i;
#line 261 "claic1.f"
	    *sestpr = sqrt(t + 1.) * absest;
#line 262 "claic1.f"
	    return 0;
#line 263 "claic1.f"
	}

#line 265 "claic1.f"
    } else if (*job == 2) {

/*        Estimating smallest singular value */

/*        special cases */

#line 271 "claic1.f"
	if (*sest == 0.) {
#line 272 "claic1.f"
	    *sestpr = 0.;
#line 273 "claic1.f"
	    if (max(absgam,absalp) == 0.) {
#line 274 "claic1.f"
		sine.r = 1., sine.i = 0.;
#line 275 "claic1.f"
		cosine.r = 0., cosine.i = 0.;
#line 276 "claic1.f"
	    } else {
#line 277 "claic1.f"
		d_cnjg(&z__2, gamma);
#line 277 "claic1.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 277 "claic1.f"
		sine.r = z__1.r, sine.i = z__1.i;
#line 278 "claic1.f"
		d_cnjg(&z__1, &alpha);
#line 278 "claic1.f"
		cosine.r = z__1.r, cosine.i = z__1.i;
#line 279 "claic1.f"
	    }
/* Computing MAX */
#line 280 "claic1.f"
	    d__1 = z_abs(&sine), d__2 = z_abs(&cosine);
#line 280 "claic1.f"
	    s1 = max(d__1,d__2);
#line 281 "claic1.f"
	    z__1.r = sine.r / s1, z__1.i = sine.i / s1;
#line 281 "claic1.f"
	    s->r = z__1.r, s->i = z__1.i;
#line 282 "claic1.f"
	    z__1.r = cosine.r / s1, z__1.i = cosine.i / s1;
#line 282 "claic1.f"
	    c__->r = z__1.r, c__->i = z__1.i;
#line 283 "claic1.f"
	    d_cnjg(&z__4, s);
#line 283 "claic1.f"
	    z__3.r = s->r * z__4.r - s->i * z__4.i, z__3.i = s->r * z__4.i + 
		    s->i * z__4.r;
#line 283 "claic1.f"
	    d_cnjg(&z__6, c__);
#line 283 "claic1.f"
	    z__5.r = c__->r * z__6.r - c__->i * z__6.i, z__5.i = c__->r * 
		    z__6.i + c__->i * z__6.r;
#line 283 "claic1.f"
	    z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
#line 283 "claic1.f"
	    z_sqrt(&z__1, &z__2);
#line 283 "claic1.f"
	    tmp = z__1.r;
#line 284 "claic1.f"
	    z__1.r = s->r / tmp, z__1.i = s->i / tmp;
#line 284 "claic1.f"
	    s->r = z__1.r, s->i = z__1.i;
#line 285 "claic1.f"
	    z__1.r = c__->r / tmp, z__1.i = c__->i / tmp;
#line 285 "claic1.f"
	    c__->r = z__1.r, c__->i = z__1.i;
#line 286 "claic1.f"
	    return 0;
#line 287 "claic1.f"
	} else if (absgam <= eps * absest) {
#line 288 "claic1.f"
	    s->r = 0., s->i = 0.;
#line 289 "claic1.f"
	    c__->r = 1., c__->i = 0.;
#line 290 "claic1.f"
	    *sestpr = absgam;
#line 291 "claic1.f"
	    return 0;
#line 292 "claic1.f"
	} else if (absalp <= eps * absest) {
#line 293 "claic1.f"
	    s1 = absgam;
#line 294 "claic1.f"
	    s2 = absest;
#line 295 "claic1.f"
	    if (s1 <= s2) {
#line 296 "claic1.f"
		s->r = 0., s->i = 0.;
#line 297 "claic1.f"
		c__->r = 1., c__->i = 0.;
#line 298 "claic1.f"
		*sestpr = s1;
#line 299 "claic1.f"
	    } else {
#line 300 "claic1.f"
		s->r = 1., s->i = 0.;
#line 301 "claic1.f"
		c__->r = 0., c__->i = 0.;
#line 302 "claic1.f"
		*sestpr = s2;
#line 303 "claic1.f"
	    }
#line 304 "claic1.f"
	    return 0;
#line 305 "claic1.f"
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
#line 306 "claic1.f"
	    s1 = absgam;
#line 307 "claic1.f"
	    s2 = absalp;
#line 308 "claic1.f"
	    if (s1 <= s2) {
#line 309 "claic1.f"
		tmp = s1 / s2;
#line 310 "claic1.f"
		scl = sqrt(tmp * tmp + 1.);
#line 311 "claic1.f"
		*sestpr = absest * (tmp / scl);
#line 312 "claic1.f"
		d_cnjg(&z__4, gamma);
#line 312 "claic1.f"
		z__3.r = z__4.r / s2, z__3.i = z__4.i / s2;
#line 312 "claic1.f"
		z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 312 "claic1.f"
		z__1.r = z__2.r / scl, z__1.i = z__2.i / scl;
#line 312 "claic1.f"
		s->r = z__1.r, s->i = z__1.i;
#line 313 "claic1.f"
		d_cnjg(&z__3, &alpha);
#line 313 "claic1.f"
		z__2.r = z__3.r / s2, z__2.i = z__3.i / s2;
#line 313 "claic1.f"
		z__1.r = z__2.r / scl, z__1.i = z__2.i / scl;
#line 313 "claic1.f"
		c__->r = z__1.r, c__->i = z__1.i;
#line 314 "claic1.f"
	    } else {
#line 315 "claic1.f"
		tmp = s2 / s1;
#line 316 "claic1.f"
		scl = sqrt(tmp * tmp + 1.);
#line 317 "claic1.f"
		*sestpr = absest / scl;
#line 318 "claic1.f"
		d_cnjg(&z__4, gamma);
#line 318 "claic1.f"
		z__3.r = z__4.r / s1, z__3.i = z__4.i / s1;
#line 318 "claic1.f"
		z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 318 "claic1.f"
		z__1.r = z__2.r / scl, z__1.i = z__2.i / scl;
#line 318 "claic1.f"
		s->r = z__1.r, s->i = z__1.i;
#line 319 "claic1.f"
		d_cnjg(&z__3, &alpha);
#line 319 "claic1.f"
		z__2.r = z__3.r / s1, z__2.i = z__3.i / s1;
#line 319 "claic1.f"
		z__1.r = z__2.r / scl, z__1.i = z__2.i / scl;
#line 319 "claic1.f"
		c__->r = z__1.r, c__->i = z__1.i;
#line 320 "claic1.f"
	    }
#line 321 "claic1.f"
	    return 0;
#line 322 "claic1.f"
	} else {

/*           normal case */

#line 326 "claic1.f"
	    zeta1 = absalp / absest;
#line 327 "claic1.f"
	    zeta2 = absgam / absest;

/* Computing MAX */
#line 329 "claic1.f"
	    d__1 = zeta1 * zeta1 + 1. + zeta1 * zeta2, d__2 = zeta1 * zeta2 + 
		    zeta2 * zeta2;
#line 329 "claic1.f"
	    norma = max(d__1,d__2);

/*           See if root is closer to zero or to ONE */

#line 334 "claic1.f"
	    test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
#line 335 "claic1.f"
	    if (test >= 0.) {

/*              root is close to zero, compute directly */

#line 339 "claic1.f"
		b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
#line 340 "claic1.f"
		d__1 = zeta2 * zeta2;
#line 340 "claic1.f"
		c__->r = d__1, c__->i = 0.;
#line 341 "claic1.f"
		d__2 = b * b;
#line 341 "claic1.f"
		z__2.r = d__2 - c__->r, z__2.i = -c__->i;
#line 341 "claic1.f"
		d__1 = b + sqrt(z_abs(&z__2));
#line 341 "claic1.f"
		z__1.r = c__->r / d__1, z__1.i = c__->i / d__1;
#line 341 "claic1.f"
		t = z__1.r;
#line 342 "claic1.f"
		z__2.r = alpha.r / absest, z__2.i = alpha.i / absest;
#line 342 "claic1.f"
		d__1 = 1. - t;
#line 342 "claic1.f"
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 342 "claic1.f"
		sine.r = z__1.r, sine.i = z__1.i;
#line 343 "claic1.f"
		z__3.r = gamma->r / absest, z__3.i = gamma->i / absest;
#line 343 "claic1.f"
		z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 343 "claic1.f"
		z__1.r = z__2.r / t, z__1.i = z__2.i / t;
#line 343 "claic1.f"
		cosine.r = z__1.r, cosine.i = z__1.i;
#line 344 "claic1.f"
		*sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
#line 345 "claic1.f"
	    } else {

/*              root is closer to ONE, shift by that amount */

#line 349 "claic1.f"
		b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
#line 350 "claic1.f"
		d__1 = zeta1 * zeta1;
#line 350 "claic1.f"
		c__->r = d__1, c__->i = 0.;
#line 351 "claic1.f"
		if (b >= 0.) {
#line 352 "claic1.f"
		    z__2.r = -c__->r, z__2.i = -c__->i;
#line 352 "claic1.f"
		    d__1 = b * b;
#line 352 "claic1.f"
		    z__5.r = d__1 + c__->r, z__5.i = c__->i;
#line 352 "claic1.f"
		    z_sqrt(&z__4, &z__5);
#line 352 "claic1.f"
		    z__3.r = b + z__4.r, z__3.i = z__4.i;
#line 352 "claic1.f"
		    z_div(&z__1, &z__2, &z__3);
#line 352 "claic1.f"
		    t = z__1.r;
#line 353 "claic1.f"
		} else {
#line 354 "claic1.f"
		    d__1 = b * b;
#line 354 "claic1.f"
		    z__3.r = d__1 + c__->r, z__3.i = c__->i;
#line 354 "claic1.f"
		    z_sqrt(&z__2, &z__3);
#line 354 "claic1.f"
		    z__1.r = b - z__2.r, z__1.i = -z__2.i;
#line 354 "claic1.f"
		    t = z__1.r;
#line 355 "claic1.f"
		}
#line 356 "claic1.f"
		z__3.r = alpha.r / absest, z__3.i = alpha.i / absest;
#line 356 "claic1.f"
		z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 356 "claic1.f"
		z__1.r = z__2.r / t, z__1.i = z__2.i / t;
#line 356 "claic1.f"
		sine.r = z__1.r, sine.i = z__1.i;
#line 357 "claic1.f"
		z__3.r = gamma->r / absest, z__3.i = gamma->i / absest;
#line 357 "claic1.f"
		z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 357 "claic1.f"
		d__1 = t + 1.;
#line 357 "claic1.f"
		z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
#line 357 "claic1.f"
		cosine.r = z__1.r, cosine.i = z__1.i;
#line 358 "claic1.f"
		*sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
#line 359 "claic1.f"
	    }
#line 360 "claic1.f"
	    d_cnjg(&z__4, &sine);
#line 360 "claic1.f"
	    z__3.r = sine.r * z__4.r - sine.i * z__4.i, z__3.i = sine.r * 
		    z__4.i + sine.i * z__4.r;
#line 360 "claic1.f"
	    d_cnjg(&z__6, &cosine);
#line 360 "claic1.f"
	    z__5.r = cosine.r * z__6.r - cosine.i * z__6.i, z__5.i = cosine.r 
		    * z__6.i + cosine.i * z__6.r;
#line 360 "claic1.f"
	    z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
#line 360 "claic1.f"
	    z_sqrt(&z__1, &z__2);
#line 360 "claic1.f"
	    tmp = z__1.r;
#line 361 "claic1.f"
	    z__1.r = sine.r / tmp, z__1.i = sine.i / tmp;
#line 361 "claic1.f"
	    s->r = z__1.r, s->i = z__1.i;
#line 362 "claic1.f"
	    z__1.r = cosine.r / tmp, z__1.i = cosine.i / tmp;
#line 362 "claic1.f"
	    c__->r = z__1.r, c__->i = z__1.i;
#line 363 "claic1.f"
	    return 0;

#line 365 "claic1.f"
	}
#line 366 "claic1.f"
    }
#line 367 "claic1.f"
    return 0;

/*     End of CLAIC1 */

} /* claic1_ */

