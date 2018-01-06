#line 1 "slasv2.f"
/* slasv2.f -- translated by f2c (version 20100827).
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

#line 1 "slasv2.f"
/* Table of constant values */

static doublereal c_b3 = 2.;
static doublereal c_b4 = 1.;

/* > \brief \b SLASV2 computes the singular value decomposition of a 2-by-2 triangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASV2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasv2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasv2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasv2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL ) */

/*       .. Scalar Arguments .. */
/*       REAL               CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASV2 computes the singular value decomposition of a 2-by-2 */
/* > triangular matrix */
/* >    [  F   G  ] */
/* >    [  0   H  ]. */
/* > On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the */
/* > smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and */
/* > right singular vectors for abs(SSMAX), giving the decomposition */
/* > */
/* >    [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ] */
/* >    [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ]. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] F */
/* > \verbatim */
/* >          F is REAL */
/* >          The (1,1) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is REAL */
/* >          The (1,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* >          H is REAL */
/* >          The (2,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMIN */
/* > \verbatim */
/* >          SSMIN is REAL */
/* >          abs(SSMIN) is the smaller singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMAX */
/* > \verbatim */
/* >          SSMAX is REAL */
/* >          abs(SSMAX) is the larger singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] SNL */
/* > \verbatim */
/* >          SNL is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] CSL */
/* > \verbatim */
/* >          CSL is REAL */
/* >          The vector (CSL, SNL) is a unit left singular vector for the */
/* >          singular value abs(SSMAX). */
/* > \endverbatim */
/* > */
/* > \param[out] SNR */
/* > \verbatim */
/* >          SNR is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] CSR */
/* > \verbatim */
/* >          CSR is REAL */
/* >          The vector (CSR, SNR) is a unit right singular vector for the */
/* >          singular value abs(SSMAX). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Any input parameter may be aliased with any output parameter. */
/* > */
/* >  Barring over/underflow and assuming a guard digit in subtraction, all */
/* >  output quantities are correct to within a few units in the last */
/* >  place (ulps). */
/* > */
/* >  In IEEE arithmetic, the code works correctly if one matrix element is */
/* >  infinite. */
/* > */
/* >  Overflow will not occur unless the largest singular value itself */
/* >  overflows or is within a few ulps of overflow. (On machines with */
/* >  partial overflow, like the Cray, overflow may occur if the largest */
/* >  singular value is within a factor of 2 of overflow.) */
/* > */
/* >  Underflow is harmless if underflow is gradual. Otherwise, results */
/* >  may correspond to a matrix modified by perturbations of size near */
/* >  the underflow threshold. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slasv2_(doublereal *f, doublereal *g, doublereal *h__, 
	doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *
	csr, doublereal *snl, doublereal *csl)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal a, d__, l, m, r__, s, t, fa, ga, ha, ft, gt, ht, mm, tt,
	     clt, crt, slt, srt;
    static integer pmax;
    static doublereal temp;
    static logical swap;
    static doublereal tsign;
    static logical gasmal;
    extern doublereal slamch_(char *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 179 "slasv2.f"
    ft = *f;
#line 180 "slasv2.f"
    fa = abs(ft);
#line 181 "slasv2.f"
    ht = *h__;
#line 182 "slasv2.f"
    ha = abs(*h__);

/*     PMAX points to the maximum absolute element of matrix */
/*       PMAX = 1 if F largest in absolute values */
/*       PMAX = 2 if G largest in absolute values */
/*       PMAX = 3 if H largest in absolute values */

#line 189 "slasv2.f"
    pmax = 1;
#line 190 "slasv2.f"
    swap = ha > fa;
#line 191 "slasv2.f"
    if (swap) {
#line 192 "slasv2.f"
	pmax = 3;
#line 193 "slasv2.f"
	temp = ft;
#line 194 "slasv2.f"
	ft = ht;
#line 195 "slasv2.f"
	ht = temp;
#line 196 "slasv2.f"
	temp = fa;
#line 197 "slasv2.f"
	fa = ha;
#line 198 "slasv2.f"
	ha = temp;

/*        Now FA .ge. HA */

#line 202 "slasv2.f"
    }
#line 203 "slasv2.f"
    gt = *g;
#line 204 "slasv2.f"
    ga = abs(gt);
#line 205 "slasv2.f"
    if (ga == 0.) {

/*        Diagonal matrix */

#line 209 "slasv2.f"
	*ssmin = ha;
#line 210 "slasv2.f"
	*ssmax = fa;
#line 211 "slasv2.f"
	clt = 1.;
#line 212 "slasv2.f"
	crt = 1.;
#line 213 "slasv2.f"
	slt = 0.;
#line 214 "slasv2.f"
	srt = 0.;
#line 215 "slasv2.f"
    } else {
#line 216 "slasv2.f"
	gasmal = TRUE_;
#line 217 "slasv2.f"
	if (ga > fa) {
#line 218 "slasv2.f"
	    pmax = 2;
#line 219 "slasv2.f"
	    if (fa / ga < slamch_("EPS", (ftnlen)3)) {

/*              Case of very large GA */

#line 223 "slasv2.f"
		gasmal = FALSE_;
#line 224 "slasv2.f"
		*ssmax = ga;
#line 225 "slasv2.f"
		if (ha > 1.) {
#line 226 "slasv2.f"
		    *ssmin = fa / (ga / ha);
#line 227 "slasv2.f"
		} else {
#line 228 "slasv2.f"
		    *ssmin = fa / ga * ha;
#line 229 "slasv2.f"
		}
#line 230 "slasv2.f"
		clt = 1.;
#line 231 "slasv2.f"
		slt = ht / gt;
#line 232 "slasv2.f"
		srt = 1.;
#line 233 "slasv2.f"
		crt = ft / gt;
#line 234 "slasv2.f"
	    }
#line 235 "slasv2.f"
	}
#line 236 "slasv2.f"
	if (gasmal) {

/*           Normal case */

#line 240 "slasv2.f"
	    d__ = fa - ha;
#line 241 "slasv2.f"
	    if (d__ == fa) {

/*              Copes with infinite F or H */

#line 245 "slasv2.f"
		l = 1.;
#line 246 "slasv2.f"
	    } else {
#line 247 "slasv2.f"
		l = d__ / fa;
#line 248 "slasv2.f"
	    }

/*           Note that 0 .le. L .le. 1 */

#line 252 "slasv2.f"
	    m = gt / ft;

/*           Note that abs(M) .le. 1/macheps */

#line 256 "slasv2.f"
	    t = 2. - l;

/*           Note that T .ge. 1 */

#line 260 "slasv2.f"
	    mm = m * m;
#line 261 "slasv2.f"
	    tt = t * t;
#line 262 "slasv2.f"
	    s = sqrt(tt + mm);

/*           Note that 1 .le. S .le. 1 + 1/macheps */

#line 266 "slasv2.f"
	    if (l == 0.) {
#line 267 "slasv2.f"
		r__ = abs(m);
#line 268 "slasv2.f"
	    } else {
#line 269 "slasv2.f"
		r__ = sqrt(l * l + mm);
#line 270 "slasv2.f"
	    }

/*           Note that 0 .le. R .le. 1 + 1/macheps */

#line 274 "slasv2.f"
	    a = (s + r__) * .5;

/*           Note that 1 .le. A .le. 1 + abs(M) */

#line 278 "slasv2.f"
	    *ssmin = ha / a;
#line 279 "slasv2.f"
	    *ssmax = fa * a;
#line 280 "slasv2.f"
	    if (mm == 0.) {

/*              Note that M is very tiny */

#line 284 "slasv2.f"
		if (l == 0.) {
#line 285 "slasv2.f"
		    t = d_sign(&c_b3, &ft) * d_sign(&c_b4, &gt);
#line 286 "slasv2.f"
		} else {
#line 287 "slasv2.f"
		    t = gt / d_sign(&d__, &ft) + m / t;
#line 288 "slasv2.f"
		}
#line 289 "slasv2.f"
	    } else {
#line 290 "slasv2.f"
		t = (m / (s + t) + m / (r__ + l)) * (a + 1.);
#line 291 "slasv2.f"
	    }
#line 292 "slasv2.f"
	    l = sqrt(t * t + 4.);
#line 293 "slasv2.f"
	    crt = 2. / l;
#line 294 "slasv2.f"
	    srt = t / l;
#line 295 "slasv2.f"
	    clt = (crt + srt * m) / a;
#line 296 "slasv2.f"
	    slt = ht / ft * srt / a;
#line 297 "slasv2.f"
	}
#line 298 "slasv2.f"
    }
#line 299 "slasv2.f"
    if (swap) {
#line 300 "slasv2.f"
	*csl = srt;
#line 301 "slasv2.f"
	*snl = crt;
#line 302 "slasv2.f"
	*csr = slt;
#line 303 "slasv2.f"
	*snr = clt;
#line 304 "slasv2.f"
    } else {
#line 305 "slasv2.f"
	*csl = clt;
#line 306 "slasv2.f"
	*snl = slt;
#line 307 "slasv2.f"
	*csr = crt;
#line 308 "slasv2.f"
	*snr = srt;
#line 309 "slasv2.f"
    }

/*     Correct signs of SSMAX and SSMIN */

#line 313 "slasv2.f"
    if (pmax == 1) {
#line 313 "slasv2.f"
	tsign = d_sign(&c_b4, csr) * d_sign(&c_b4, csl) * d_sign(&c_b4, f);
#line 313 "slasv2.f"
    }
#line 315 "slasv2.f"
    if (pmax == 2) {
#line 315 "slasv2.f"
	tsign = d_sign(&c_b4, snr) * d_sign(&c_b4, csl) * d_sign(&c_b4, g);
#line 315 "slasv2.f"
    }
#line 317 "slasv2.f"
    if (pmax == 3) {
#line 317 "slasv2.f"
	tsign = d_sign(&c_b4, snr) * d_sign(&c_b4, snl) * d_sign(&c_b4, h__);
#line 317 "slasv2.f"
    }
#line 319 "slasv2.f"
    *ssmax = d_sign(ssmax, &tsign);
#line 320 "slasv2.f"
    d__1 = tsign * d_sign(&c_b4, f) * d_sign(&c_b4, h__);
#line 320 "slasv2.f"
    *ssmin = d_sign(ssmin, &d__1);
#line 321 "slasv2.f"
    return 0;

/*     End of SLASV2 */

} /* slasv2_ */

