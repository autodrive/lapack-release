#line 1 "slanv2.f"
/* slanv2.f -- translated by f2c (version 20100827).
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

#line 1 "slanv2.f"
/* Table of constant values */

static doublereal c_b4 = 1.;

/* > \brief \b SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric matrix in standard form. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANV2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slanv2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slanv2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slanv2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN ) */

/*       .. Scalar Arguments .. */
/*       REAL               A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric */
/* > matrix in standard form: */
/* > */
/* >      [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ] */
/* >      [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ] */
/* > */
/* > where either */
/* > 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or */
/* > 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex */
/* > conjugate eigenvalues. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL */
/* >          On entry, the elements of the input matrix. */
/* >          On exit, they are overwritten by the elements of the */
/* >          standardised Schur form. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1R */
/* > \verbatim */
/* >          RT1R is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] RT1I */
/* > \verbatim */
/* >          RT1I is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] RT2R */
/* > \verbatim */
/* >          RT2R is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] RT2I */
/* > \verbatim */
/* >          RT2I is REAL */
/* >          The real and imaginary parts of the eigenvalues. If the */
/* >          eigenvalues are a complex conjugate pair, RT1I > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* >          CS is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* >          SN is REAL */
/* >          Parameters of the rotation matrix. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified by V. Sima, Research Institute for Informatics, Bucharest, */
/* >  Romania, to reduce the risk of cancellation errors, */
/* >  when computing real eigenvalues, and to ensure, if possible, that */
/* >  abs(RT1R) >= abs(RT2R). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slanv2_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r,
	 doublereal *rt2i, doublereal *cs, doublereal *sn)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal p, z__, aa, bb, cc, dd, cs1, sn1, sab, sac, eps, tau, 
	    temp, scale, bcmax, bcmis, sigma;
    extern doublereal slapy2_(doublereal *, doublereal *), slamch_(char *, 
	    ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 160 "slanv2.f"
    eps = slamch_("P", (ftnlen)1);
#line 161 "slanv2.f"
    if (*c__ == 0.) {
#line 162 "slanv2.f"
	*cs = 1.;
#line 163 "slanv2.f"
	*sn = 0.;
#line 164 "slanv2.f"
	goto L10;

#line 166 "slanv2.f"
    } else if (*b == 0.) {

/*        Swap rows and columns */

#line 170 "slanv2.f"
	*cs = 0.;
#line 171 "slanv2.f"
	*sn = 1.;
#line 172 "slanv2.f"
	temp = *d__;
#line 173 "slanv2.f"
	*d__ = *a;
#line 174 "slanv2.f"
	*a = temp;
#line 175 "slanv2.f"
	*b = -(*c__);
#line 176 "slanv2.f"
	*c__ = 0.;
#line 177 "slanv2.f"
	goto L10;
#line 178 "slanv2.f"
    } else if (*a - *d__ == 0. && d_sign(&c_b4, b) != d_sign(&c_b4, c__)) {
#line 180 "slanv2.f"
	*cs = 1.;
#line 181 "slanv2.f"
	*sn = 0.;
#line 182 "slanv2.f"
	goto L10;
#line 183 "slanv2.f"
    } else {

#line 185 "slanv2.f"
	temp = *a - *d__;
#line 186 "slanv2.f"
	p = temp * .5;
/* Computing MAX */
#line 187 "slanv2.f"
	d__1 = abs(*b), d__2 = abs(*c__);
#line 187 "slanv2.f"
	bcmax = max(d__1,d__2);
/* Computing MIN */
#line 188 "slanv2.f"
	d__1 = abs(*b), d__2 = abs(*c__);
#line 188 "slanv2.f"
	bcmis = min(d__1,d__2) * d_sign(&c_b4, b) * d_sign(&c_b4, c__);
/* Computing MAX */
#line 189 "slanv2.f"
	d__1 = abs(p);
#line 189 "slanv2.f"
	scale = max(d__1,bcmax);
#line 190 "slanv2.f"
	z__ = p / scale * p + bcmax / scale * bcmis;

/*        If Z is of the order of the machine accuracy, postpone the */
/*        decision on the nature of eigenvalues */

#line 195 "slanv2.f"
	if (z__ >= eps * 4.) {

/*           Real eigenvalues. Compute A and D. */

#line 199 "slanv2.f"
	    d__1 = sqrt(scale) * sqrt(z__);
#line 199 "slanv2.f"
	    z__ = p + d_sign(&d__1, &p);
#line 200 "slanv2.f"
	    *a = *d__ + z__;
#line 201 "slanv2.f"
	    *d__ -= bcmax / z__ * bcmis;

/*           Compute B and the rotation matrix */

#line 205 "slanv2.f"
	    tau = slapy2_(c__, &z__);
#line 206 "slanv2.f"
	    *cs = z__ / tau;
#line 207 "slanv2.f"
	    *sn = *c__ / tau;
#line 208 "slanv2.f"
	    *b -= *c__;
#line 209 "slanv2.f"
	    *c__ = 0.;
#line 210 "slanv2.f"
	} else {

/*           Complex eigenvalues, or real (almost) equal eigenvalues. */
/*           Make diagonal elements equal. */

#line 215 "slanv2.f"
	    sigma = *b + *c__;
#line 216 "slanv2.f"
	    tau = slapy2_(&sigma, &temp);
#line 217 "slanv2.f"
	    *cs = sqrt((abs(sigma) / tau + 1.) * .5);
#line 218 "slanv2.f"
	    *sn = -(p / (tau * *cs)) * d_sign(&c_b4, &sigma);

/*           Compute [ AA  BB ] = [ A  B ] [ CS -SN ] */
/*                   [ CC  DD ]   [ C  D ] [ SN  CS ] */

#line 223 "slanv2.f"
	    aa = *a * *cs + *b * *sn;
#line 224 "slanv2.f"
	    bb = -(*a) * *sn + *b * *cs;
#line 225 "slanv2.f"
	    cc = *c__ * *cs + *d__ * *sn;
#line 226 "slanv2.f"
	    dd = -(*c__) * *sn + *d__ * *cs;

/*           Compute [ A  B ] = [ CS  SN ] [ AA  BB ] */
/*                   [ C  D ]   [-SN  CS ] [ CC  DD ] */

#line 231 "slanv2.f"
	    *a = aa * *cs + cc * *sn;
#line 232 "slanv2.f"
	    *b = bb * *cs + dd * *sn;
#line 233 "slanv2.f"
	    *c__ = -aa * *sn + cc * *cs;
#line 234 "slanv2.f"
	    *d__ = -bb * *sn + dd * *cs;

#line 236 "slanv2.f"
	    temp = (*a + *d__) * .5;
#line 237 "slanv2.f"
	    *a = temp;
#line 238 "slanv2.f"
	    *d__ = temp;

#line 240 "slanv2.f"
	    if (*c__ != 0.) {
#line 241 "slanv2.f"
		if (*b != 0.) {
#line 242 "slanv2.f"
		    if (d_sign(&c_b4, b) == d_sign(&c_b4, c__)) {

/*                    Real eigenvalues: reduce to upper triangular form */

#line 246 "slanv2.f"
			sab = sqrt((abs(*b)));
#line 247 "slanv2.f"
			sac = sqrt((abs(*c__)));
#line 248 "slanv2.f"
			d__1 = sab * sac;
#line 248 "slanv2.f"
			p = d_sign(&d__1, c__);
#line 249 "slanv2.f"
			tau = 1. / sqrt((d__1 = *b + *c__, abs(d__1)));
#line 250 "slanv2.f"
			*a = temp + p;
#line 251 "slanv2.f"
			*d__ = temp - p;
#line 252 "slanv2.f"
			*b -= *c__;
#line 253 "slanv2.f"
			*c__ = 0.;
#line 254 "slanv2.f"
			cs1 = sab * tau;
#line 255 "slanv2.f"
			sn1 = sac * tau;
#line 256 "slanv2.f"
			temp = *cs * cs1 - *sn * sn1;
#line 257 "slanv2.f"
			*sn = *cs * sn1 + *sn * cs1;
#line 258 "slanv2.f"
			*cs = temp;
#line 259 "slanv2.f"
		    }
#line 260 "slanv2.f"
		} else {
#line 261 "slanv2.f"
		    *b = -(*c__);
#line 262 "slanv2.f"
		    *c__ = 0.;
#line 263 "slanv2.f"
		    temp = *cs;
#line 264 "slanv2.f"
		    *cs = -(*sn);
#line 265 "slanv2.f"
		    *sn = temp;
#line 266 "slanv2.f"
		}
#line 267 "slanv2.f"
	    }
#line 268 "slanv2.f"
	}

#line 270 "slanv2.f"
    }

#line 272 "slanv2.f"
L10:

/*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */

#line 276 "slanv2.f"
    *rt1r = *a;
#line 277 "slanv2.f"
    *rt2r = *d__;
#line 278 "slanv2.f"
    if (*c__ == 0.) {
#line 279 "slanv2.f"
	*rt1i = 0.;
#line 280 "slanv2.f"
	*rt2i = 0.;
#line 281 "slanv2.f"
    } else {
#line 282 "slanv2.f"
	*rt1i = sqrt((abs(*b))) * sqrt((abs(*c__)));
#line 283 "slanv2.f"
	*rt2i = -(*rt1i);
#line 284 "slanv2.f"
    }
#line 285 "slanv2.f"
    return 0;

/*     End of SLANV2 */

} /* slanv2_ */

