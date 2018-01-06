#line 1 "zlartg.f"
/* zlartg.f -- translated by f2c (version 20100827).
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

#line 1 "zlartg.f"
/* > \brief \b ZLARTG generates a plane rotation with real cosine and complex sine. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLARTG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlartg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlartg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlartg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARTG( F, G, CS, SN, R ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   CS */
/*       COMPLEX*16         F, G, R, SN */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARTG generates a plane rotation so that */
/* > */
/* >    [  CS  SN  ]     [ F ]     [ R ] */
/* >    [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1. */
/* >    [ -SN  CS  ]     [ G ]     [ 0 ] */
/* > */
/* > This is a faster version of the BLAS1 routine ZROTG, except for */
/* > the following differences: */
/* >    F and G are unchanged on return. */
/* >    If G=0, then CS=1 and SN=0. */
/* >    If F=0, then CS=0 and SN is chosen so that R is real. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] F */
/* > \verbatim */
/* >          F is COMPLEX*16 */
/* >          The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is COMPLEX*16 */
/* >          The second component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* >          CS is DOUBLE PRECISION */
/* >          The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* >          SN is COMPLEX*16 */
/* >          The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* >          R is COMPLEX*16 */
/* >          The nonzero component of the rotated vector. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel */
/* > */
/* >  This version has a few statements commented out for thread safety */
/* >  (machine parameters are computed on each entry). 10 feb 03, SJH. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlartg_(doublecomplex *f, doublecomplex *g, doublereal *
	cs, doublecomplex *sn, doublecomplex *r__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), d_imag(
	    doublecomplex *), z_abs(doublecomplex *), sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal d__;
    static integer i__;
    static doublereal f2, g2;
    static doublecomplex ff;
    static doublereal di, dr;
    static doublecomplex fs, gs;
    static doublereal f2s, g2s, eps, scale;
    static integer count;
    static doublereal safmn2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal safmx2;
    extern doublereal dlamch_(char *, ftnlen);
    extern logical disnan_(doublereal *);
    static doublereal safmin;


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
/*     LOGICAL            FIRST */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 149 "zlartg.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 150 "zlartg.f"
    eps = dlamch_("E", (ftnlen)1);
#line 151 "zlartg.f"
    d__1 = dlamch_("B", (ftnlen)1);
#line 151 "zlartg.f"
    i__1 = (integer) (log(safmin / eps) / log(dlamch_("B", (ftnlen)1)) / 2.);
#line 151 "zlartg.f"
    safmn2 = pow_di(&d__1, &i__1);
#line 153 "zlartg.f"
    safmx2 = 1. / safmn2;
/* Computing MAX */
/* Computing MAX */
#line 154 "zlartg.f"
    d__7 = (d__1 = f->r, abs(d__1)), d__8 = (d__2 = d_imag(f), abs(d__2));
/* Computing MAX */
#line 154 "zlartg.f"
    d__9 = (d__3 = g->r, abs(d__3)), d__10 = (d__4 = d_imag(g), abs(d__4));
#line 154 "zlartg.f"
    d__5 = max(d__7,d__8), d__6 = max(d__9,d__10);
#line 154 "zlartg.f"
    scale = max(d__5,d__6);
#line 155 "zlartg.f"
    fs.r = f->r, fs.i = f->i;
#line 156 "zlartg.f"
    gs.r = g->r, gs.i = g->i;
#line 157 "zlartg.f"
    count = 0;
#line 158 "zlartg.f"
    if (scale >= safmx2) {
#line 159 "zlartg.f"
L10:
#line 160 "zlartg.f"
	++count;
#line 161 "zlartg.f"
	z__1.r = safmn2 * fs.r, z__1.i = safmn2 * fs.i;
#line 161 "zlartg.f"
	fs.r = z__1.r, fs.i = z__1.i;
#line 162 "zlartg.f"
	z__1.r = safmn2 * gs.r, z__1.i = safmn2 * gs.i;
#line 162 "zlartg.f"
	gs.r = z__1.r, gs.i = z__1.i;
#line 163 "zlartg.f"
	scale *= safmn2;
#line 164 "zlartg.f"
	if (scale >= safmx2) {
#line 164 "zlartg.f"
	    goto L10;
#line 164 "zlartg.f"
	}
#line 166 "zlartg.f"
    } else if (scale <= safmn2) {
#line 167 "zlartg.f"
	d__1 = z_abs(g);
#line 167 "zlartg.f"
	if (g->r == 0. && g->i == 0. || disnan_(&d__1)) {
#line 168 "zlartg.f"
	    *cs = 1.;
#line 169 "zlartg.f"
	    sn->r = 0., sn->i = 0.;
#line 170 "zlartg.f"
	    r__->r = f->r, r__->i = f->i;
#line 171 "zlartg.f"
	    return 0;
#line 172 "zlartg.f"
	}
#line 173 "zlartg.f"
L20:
#line 174 "zlartg.f"
	--count;
#line 175 "zlartg.f"
	z__1.r = safmx2 * fs.r, z__1.i = safmx2 * fs.i;
#line 175 "zlartg.f"
	fs.r = z__1.r, fs.i = z__1.i;
#line 176 "zlartg.f"
	z__1.r = safmx2 * gs.r, z__1.i = safmx2 * gs.i;
#line 176 "zlartg.f"
	gs.r = z__1.r, gs.i = z__1.i;
#line 177 "zlartg.f"
	scale *= safmx2;
#line 178 "zlartg.f"
	if (scale <= safmn2) {
#line 178 "zlartg.f"
	    goto L20;
#line 178 "zlartg.f"
	}
#line 180 "zlartg.f"
    }
/* Computing 2nd power */
#line 181 "zlartg.f"
    d__1 = fs.r;
/* Computing 2nd power */
#line 181 "zlartg.f"
    d__2 = d_imag(&fs);
#line 181 "zlartg.f"
    f2 = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
#line 182 "zlartg.f"
    d__1 = gs.r;
/* Computing 2nd power */
#line 182 "zlartg.f"
    d__2 = d_imag(&gs);
#line 182 "zlartg.f"
    g2 = d__1 * d__1 + d__2 * d__2;
#line 183 "zlartg.f"
    if (f2 <= max(g2,1.) * safmin) {

/*        This is a rare case: F is very small. */

#line 187 "zlartg.f"
	if (f->r == 0. && f->i == 0.) {
#line 188 "zlartg.f"
	    *cs = 0.;
#line 189 "zlartg.f"
	    d__2 = g->r;
#line 189 "zlartg.f"
	    d__3 = d_imag(g);
#line 189 "zlartg.f"
	    d__1 = dlapy2_(&d__2, &d__3);
#line 189 "zlartg.f"
	    r__->r = d__1, r__->i = 0.;
/*           Do complex/real division explicitly with two real divisions */
#line 191 "zlartg.f"
	    d__1 = gs.r;
#line 191 "zlartg.f"
	    d__2 = d_imag(&gs);
#line 191 "zlartg.f"
	    d__ = dlapy2_(&d__1, &d__2);
#line 192 "zlartg.f"
	    d__1 = gs.r / d__;
#line 192 "zlartg.f"
	    d__2 = -d_imag(&gs) / d__;
#line 192 "zlartg.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 192 "zlartg.f"
	    sn->r = z__1.r, sn->i = z__1.i;
#line 193 "zlartg.f"
	    return 0;
#line 194 "zlartg.f"
	}
#line 195 "zlartg.f"
	d__1 = fs.r;
#line 195 "zlartg.f"
	d__2 = d_imag(&fs);
#line 195 "zlartg.f"
	f2s = dlapy2_(&d__1, &d__2);
/*        G2 and G2S are accurate */
/*        G2 is at least SAFMIN, and G2S is at least SAFMN2 */
#line 198 "zlartg.f"
	g2s = sqrt(g2);
/*        Error in CS from underflow in F2S is at most */
/*        UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS */
/*        If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN, */
/*        and so CS .lt. sqrt(SAFMIN) */
/*        If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN */
/*        and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS) */
/*        Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S */
#line 206 "zlartg.f"
	*cs = f2s / g2s;
/*        Make sure abs(FF) = 1 */
/*        Do complex/real division explicitly with 2 real divisions */
/* Computing MAX */
#line 209 "zlartg.f"
	d__3 = (d__1 = f->r, abs(d__1)), d__4 = (d__2 = d_imag(f), abs(d__2));
#line 209 "zlartg.f"
	if (max(d__3,d__4) > 1.) {
#line 210 "zlartg.f"
	    d__1 = f->r;
#line 210 "zlartg.f"
	    d__2 = d_imag(f);
#line 210 "zlartg.f"
	    d__ = dlapy2_(&d__1, &d__2);
#line 211 "zlartg.f"
	    d__1 = f->r / d__;
#line 211 "zlartg.f"
	    d__2 = d_imag(f) / d__;
#line 211 "zlartg.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 211 "zlartg.f"
	    ff.r = z__1.r, ff.i = z__1.i;
#line 212 "zlartg.f"
	} else {
#line 213 "zlartg.f"
	    dr = safmx2 * f->r;
#line 214 "zlartg.f"
	    di = safmx2 * d_imag(f);
#line 215 "zlartg.f"
	    d__ = dlapy2_(&dr, &di);
#line 216 "zlartg.f"
	    d__1 = dr / d__;
#line 216 "zlartg.f"
	    d__2 = di / d__;
#line 216 "zlartg.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 216 "zlartg.f"
	    ff.r = z__1.r, ff.i = z__1.i;
#line 217 "zlartg.f"
	}
#line 218 "zlartg.f"
	d__1 = gs.r / g2s;
#line 218 "zlartg.f"
	d__2 = -d_imag(&gs) / g2s;
#line 218 "zlartg.f"
	z__2.r = d__1, z__2.i = d__2;
#line 218 "zlartg.f"
	z__1.r = ff.r * z__2.r - ff.i * z__2.i, z__1.i = ff.r * z__2.i + ff.i 
		* z__2.r;
#line 218 "zlartg.f"
	sn->r = z__1.r, sn->i = z__1.i;
#line 219 "zlartg.f"
	z__2.r = *cs * f->r, z__2.i = *cs * f->i;
#line 219 "zlartg.f"
	z__3.r = sn->r * g->r - sn->i * g->i, z__3.i = sn->r * g->i + sn->i * 
		g->r;
#line 219 "zlartg.f"
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 219 "zlartg.f"
	r__->r = z__1.r, r__->i = z__1.i;
#line 220 "zlartg.f"
    } else {

/*        This is the most common case. */
/*        Neither F2 nor F2/G2 are less than SAFMIN */
/*        F2S cannot overflow, and it is accurate */

#line 226 "zlartg.f"
	f2s = sqrt(g2 / f2 + 1.);
/*        Do the F2S(real)*FS(complex) multiply with two real multiplies */
#line 228 "zlartg.f"
	d__1 = f2s * fs.r;
#line 228 "zlartg.f"
	d__2 = f2s * d_imag(&fs);
#line 228 "zlartg.f"
	z__1.r = d__1, z__1.i = d__2;
#line 228 "zlartg.f"
	r__->r = z__1.r, r__->i = z__1.i;
#line 229 "zlartg.f"
	*cs = 1. / f2s;
#line 230 "zlartg.f"
	d__ = f2 + g2;
/*        Do complex/real division explicitly with two real divisions */
#line 232 "zlartg.f"
	d__1 = r__->r / d__;
#line 232 "zlartg.f"
	d__2 = d_imag(r__) / d__;
#line 232 "zlartg.f"
	z__1.r = d__1, z__1.i = d__2;
#line 232 "zlartg.f"
	sn->r = z__1.r, sn->i = z__1.i;
#line 233 "zlartg.f"
	d_cnjg(&z__2, &gs);
#line 233 "zlartg.f"
	z__1.r = sn->r * z__2.r - sn->i * z__2.i, z__1.i = sn->r * z__2.i + 
		sn->i * z__2.r;
#line 233 "zlartg.f"
	sn->r = z__1.r, sn->i = z__1.i;
#line 234 "zlartg.f"
	if (count != 0) {
#line 235 "zlartg.f"
	    if (count > 0) {
#line 236 "zlartg.f"
		i__1 = count;
#line 236 "zlartg.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 237 "zlartg.f"
		    z__1.r = safmx2 * r__->r, z__1.i = safmx2 * r__->i;
#line 237 "zlartg.f"
		    r__->r = z__1.r, r__->i = z__1.i;
#line 238 "zlartg.f"
/* L30: */
#line 238 "zlartg.f"
		}
#line 239 "zlartg.f"
	    } else {
#line 240 "zlartg.f"
		i__1 = -count;
#line 240 "zlartg.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 241 "zlartg.f"
		    z__1.r = safmn2 * r__->r, z__1.i = safmn2 * r__->i;
#line 241 "zlartg.f"
		    r__->r = z__1.r, r__->i = z__1.i;
#line 242 "zlartg.f"
/* L40: */
#line 242 "zlartg.f"
		}
#line 243 "zlartg.f"
	    }
#line 244 "zlartg.f"
	}
#line 245 "zlartg.f"
    }
#line 246 "zlartg.f"
    return 0;

/*     End of ZLARTG */

} /* zlartg_ */

