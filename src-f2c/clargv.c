#line 1 "clargv.f"
/* clargv.f -- translated by f2c (version 20100827).
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

#line 1 "clargv.f"
/* > \brief \b CLARGV generates a vector of plane rotations with real cosines and complex sines. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clargv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clargv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clargv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARGV( N, X, INCX, Y, INCY, C, INCC ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCC, INCX, INCY, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( * ) */
/*       COMPLEX            X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARGV generates a vector of complex plane rotations with real */
/* > cosines, determined by elements of the complex vectors x and y. */
/* > For i = 1,2,...,n */
/* > */
/* >    (        c(i)   s(i) ) ( x(i) ) = ( r(i) ) */
/* >    ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  ) */
/* > */
/* >    where c(i)**2 + ABS(s(i))**2 = 1 */
/* > */
/* > The following conventions are used (these are the same as in CLARTG, */
/* > but differ from the BLAS1 routine CROTG): */
/* >    If y(i)=0, then c(i)=1 and s(i)=0. */
/* >    If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of plane rotations to be generated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (1+(N-1)*INCX) */
/* >          On entry, the vector x. */
/* >          On exit, x(i) is overwritten by r(i), for i = 1,...,n. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array, dimension (1+(N-1)*INCY) */
/* >          On entry, the vector y. */
/* >          On exit, the sines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >          The increment between elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (1+(N-1)*INCC) */
/* >          The cosines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCC */
/* > \verbatim */
/* >          INCC is INTEGER */
/* >          The increment between elements of C. INCC > 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  6-6-96 - Modified with a new algorithm by W. Kahan and J. Demmel */
/* > */
/* >  This version has a few statements commented out for thread safety */
/* >  (machine parameters are computed on each entry). 10 feb 03, SJH. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clargv_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c__, integer *incc)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), d_imag(
	    doublecomplex *), sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal d__;
    static doublecomplex f, g;
    static integer i__, j;
    static doublecomplex r__;
    static doublereal f2, g2;
    static integer ic;
    static doublereal di;
    static doublecomplex ff;
    static doublereal cs, dr;
    static doublecomplex fs, gs;
    static integer ix, iy;
    static doublecomplex sn;
    static doublereal f2s, g2s, eps, scale;
    static integer count;
    static doublereal safmn2, safmx2;
    extern doublereal slapy2_(doublereal *, doublereal *), slamch_(char *, 
	    ftnlen);
    static doublereal safmin;


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
/*     LOGICAL            FIRST */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2 */
/*     .. */
/*     .. Data statements .. */
/*     DATA               FIRST / .TRUE. / */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     IF( FIRST ) THEN */
/*        FIRST = .FALSE. */
#line 178 "clargv.f"
    /* Parameter adjustments */
#line 178 "clargv.f"
    --c__;
#line 178 "clargv.f"
    --y;
#line 178 "clargv.f"
    --x;
#line 178 "clargv.f"

#line 178 "clargv.f"
    /* Function Body */
#line 178 "clargv.f"
    safmin = slamch_("S", (ftnlen)1);
#line 179 "clargv.f"
    eps = slamch_("E", (ftnlen)1);
#line 180 "clargv.f"
    d__1 = slamch_("B", (ftnlen)1);
#line 180 "clargv.f"
    i__1 = (integer) (log(safmin / eps) / log(slamch_("B", (ftnlen)1)) / 2.);
#line 180 "clargv.f"
    safmn2 = pow_di(&d__1, &i__1);
#line 182 "clargv.f"
    safmx2 = 1. / safmn2;
/*     END IF */
#line 184 "clargv.f"
    ix = 1;
#line 185 "clargv.f"
    iy = 1;
#line 186 "clargv.f"
    ic = 1;
#line 187 "clargv.f"
    i__1 = *n;
#line 187 "clargv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 188 "clargv.f"
	i__2 = ix;
#line 188 "clargv.f"
	f.r = x[i__2].r, f.i = x[i__2].i;
#line 189 "clargv.f"
	i__2 = iy;
#line 189 "clargv.f"
	g.r = y[i__2].r, g.i = y[i__2].i;

/*        Use identical algorithm as in CLARTG */

/* Computing MAX */
/* Computing MAX */
#line 193 "clargv.f"
	d__7 = (d__1 = f.r, abs(d__1)), d__8 = (d__2 = d_imag(&f), abs(d__2));
/* Computing MAX */
#line 193 "clargv.f"
	d__9 = (d__3 = g.r, abs(d__3)), d__10 = (d__4 = d_imag(&g), abs(d__4))
		;
#line 193 "clargv.f"
	d__5 = max(d__7,d__8), d__6 = max(d__9,d__10);
#line 193 "clargv.f"
	scale = max(d__5,d__6);
#line 194 "clargv.f"
	fs.r = f.r, fs.i = f.i;
#line 195 "clargv.f"
	gs.r = g.r, gs.i = g.i;
#line 196 "clargv.f"
	count = 0;
#line 197 "clargv.f"
	if (scale >= safmx2) {
#line 198 "clargv.f"
L10:
#line 199 "clargv.f"
	    ++count;
#line 200 "clargv.f"
	    z__1.r = safmn2 * fs.r, z__1.i = safmn2 * fs.i;
#line 200 "clargv.f"
	    fs.r = z__1.r, fs.i = z__1.i;
#line 201 "clargv.f"
	    z__1.r = safmn2 * gs.r, z__1.i = safmn2 * gs.i;
#line 201 "clargv.f"
	    gs.r = z__1.r, gs.i = z__1.i;
#line 202 "clargv.f"
	    scale *= safmn2;
#line 203 "clargv.f"
	    if (scale >= safmx2) {
#line 203 "clargv.f"
		goto L10;
#line 203 "clargv.f"
	    }
#line 205 "clargv.f"
	} else if (scale <= safmn2) {
#line 206 "clargv.f"
	    if (g.r == 0. && g.i == 0.) {
#line 207 "clargv.f"
		cs = 1.;
#line 208 "clargv.f"
		sn.r = 0., sn.i = 0.;
#line 209 "clargv.f"
		r__.r = f.r, r__.i = f.i;
#line 210 "clargv.f"
		goto L50;
#line 211 "clargv.f"
	    }
#line 212 "clargv.f"
L20:
#line 213 "clargv.f"
	    --count;
#line 214 "clargv.f"
	    z__1.r = safmx2 * fs.r, z__1.i = safmx2 * fs.i;
#line 214 "clargv.f"
	    fs.r = z__1.r, fs.i = z__1.i;
#line 215 "clargv.f"
	    z__1.r = safmx2 * gs.r, z__1.i = safmx2 * gs.i;
#line 215 "clargv.f"
	    gs.r = z__1.r, gs.i = z__1.i;
#line 216 "clargv.f"
	    scale *= safmx2;
#line 217 "clargv.f"
	    if (scale <= safmn2) {
#line 217 "clargv.f"
		goto L20;
#line 217 "clargv.f"
	    }
#line 219 "clargv.f"
	}
/* Computing 2nd power */
#line 220 "clargv.f"
	d__1 = fs.r;
/* Computing 2nd power */
#line 220 "clargv.f"
	d__2 = d_imag(&fs);
#line 220 "clargv.f"
	f2 = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
#line 221 "clargv.f"
	d__1 = gs.r;
/* Computing 2nd power */
#line 221 "clargv.f"
	d__2 = d_imag(&gs);
#line 221 "clargv.f"
	g2 = d__1 * d__1 + d__2 * d__2;
#line 222 "clargv.f"
	if (f2 <= max(g2,1.) * safmin) {

/*           This is a rare case: F is very small. */

#line 226 "clargv.f"
	    if (f.r == 0. && f.i == 0.) {
#line 227 "clargv.f"
		cs = 0.;
#line 228 "clargv.f"
		d__2 = g.r;
#line 228 "clargv.f"
		d__3 = d_imag(&g);
#line 228 "clargv.f"
		d__1 = slapy2_(&d__2, &d__3);
#line 228 "clargv.f"
		r__.r = d__1, r__.i = 0.;
/*              Do complex/real division explicitly with two real */
/*              divisions */
#line 231 "clargv.f"
		d__1 = gs.r;
#line 231 "clargv.f"
		d__2 = d_imag(&gs);
#line 231 "clargv.f"
		d__ = slapy2_(&d__1, &d__2);
#line 232 "clargv.f"
		d__1 = gs.r / d__;
#line 232 "clargv.f"
		d__2 = -d_imag(&gs) / d__;
#line 232 "clargv.f"
		z__1.r = d__1, z__1.i = d__2;
#line 232 "clargv.f"
		sn.r = z__1.r, sn.i = z__1.i;
#line 233 "clargv.f"
		goto L50;
#line 234 "clargv.f"
	    }
#line 235 "clargv.f"
	    d__1 = fs.r;
#line 235 "clargv.f"
	    d__2 = d_imag(&fs);
#line 235 "clargv.f"
	    f2s = slapy2_(&d__1, &d__2);
/*           G2 and G2S are accurate */
/*           G2 is at least SAFMIN, and G2S is at least SAFMN2 */
#line 238 "clargv.f"
	    g2s = sqrt(g2);
/*           Error in CS from underflow in F2S is at most */
/*           UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS */
/*           If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN, */
/*           and so CS .lt. sqrt(SAFMIN) */
/*           If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN */
/*           and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS) */
/*           Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S */
#line 246 "clargv.f"
	    cs = f2s / g2s;
/*           Make sure abs(FF) = 1 */
/*           Do complex/real division explicitly with 2 real divisions */
/* Computing MAX */
#line 249 "clargv.f"
	    d__3 = (d__1 = f.r, abs(d__1)), d__4 = (d__2 = d_imag(&f), abs(
		    d__2));
#line 249 "clargv.f"
	    if (max(d__3,d__4) > 1.) {
#line 250 "clargv.f"
		d__1 = f.r;
#line 250 "clargv.f"
		d__2 = d_imag(&f);
#line 250 "clargv.f"
		d__ = slapy2_(&d__1, &d__2);
#line 251 "clargv.f"
		d__1 = f.r / d__;
#line 251 "clargv.f"
		d__2 = d_imag(&f) / d__;
#line 251 "clargv.f"
		z__1.r = d__1, z__1.i = d__2;
#line 251 "clargv.f"
		ff.r = z__1.r, ff.i = z__1.i;
#line 252 "clargv.f"
	    } else {
#line 253 "clargv.f"
		dr = safmx2 * f.r;
#line 254 "clargv.f"
		di = safmx2 * d_imag(&f);
#line 255 "clargv.f"
		d__ = slapy2_(&dr, &di);
#line 256 "clargv.f"
		d__1 = dr / d__;
#line 256 "clargv.f"
		d__2 = di / d__;
#line 256 "clargv.f"
		z__1.r = d__1, z__1.i = d__2;
#line 256 "clargv.f"
		ff.r = z__1.r, ff.i = z__1.i;
#line 257 "clargv.f"
	    }
#line 258 "clargv.f"
	    d__1 = gs.r / g2s;
#line 258 "clargv.f"
	    d__2 = -d_imag(&gs) / g2s;
#line 258 "clargv.f"
	    z__2.r = d__1, z__2.i = d__2;
#line 258 "clargv.f"
	    z__1.r = ff.r * z__2.r - ff.i * z__2.i, z__1.i = ff.r * z__2.i + 
		    ff.i * z__2.r;
#line 258 "clargv.f"
	    sn.r = z__1.r, sn.i = z__1.i;
#line 259 "clargv.f"
	    z__2.r = cs * f.r, z__2.i = cs * f.i;
#line 259 "clargv.f"
	    z__3.r = sn.r * g.r - sn.i * g.i, z__3.i = sn.r * g.i + sn.i * 
		    g.r;
#line 259 "clargv.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 259 "clargv.f"
	    r__.r = z__1.r, r__.i = z__1.i;
#line 260 "clargv.f"
	} else {

/*           This is the most common case. */
/*           Neither F2 nor F2/G2 are less than SAFMIN */
/*           F2S cannot overflow, and it is accurate */

#line 266 "clargv.f"
	    f2s = sqrt(g2 / f2 + 1.);
/*           Do the F2S(real)*FS(complex) multiply with two real */
/*           multiplies */
#line 269 "clargv.f"
	    d__1 = f2s * fs.r;
#line 269 "clargv.f"
	    d__2 = f2s * d_imag(&fs);
#line 269 "clargv.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 269 "clargv.f"
	    r__.r = z__1.r, r__.i = z__1.i;
#line 270 "clargv.f"
	    cs = 1. / f2s;
#line 271 "clargv.f"
	    d__ = f2 + g2;
/*           Do complex/real division explicitly with two real divisions */
#line 273 "clargv.f"
	    d__1 = r__.r / d__;
#line 273 "clargv.f"
	    d__2 = d_imag(&r__) / d__;
#line 273 "clargv.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 273 "clargv.f"
	    sn.r = z__1.r, sn.i = z__1.i;
#line 274 "clargv.f"
	    d_cnjg(&z__2, &gs);
#line 274 "clargv.f"
	    z__1.r = sn.r * z__2.r - sn.i * z__2.i, z__1.i = sn.r * z__2.i + 
		    sn.i * z__2.r;
#line 274 "clargv.f"
	    sn.r = z__1.r, sn.i = z__1.i;
#line 275 "clargv.f"
	    if (count != 0) {
#line 276 "clargv.f"
		if (count > 0) {
#line 277 "clargv.f"
		    i__2 = count;
#line 277 "clargv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 278 "clargv.f"
			z__1.r = safmx2 * r__.r, z__1.i = safmx2 * r__.i;
#line 278 "clargv.f"
			r__.r = z__1.r, r__.i = z__1.i;
#line 279 "clargv.f"
/* L30: */
#line 279 "clargv.f"
		    }
#line 280 "clargv.f"
		} else {
#line 281 "clargv.f"
		    i__2 = -count;
#line 281 "clargv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 282 "clargv.f"
			z__1.r = safmn2 * r__.r, z__1.i = safmn2 * r__.i;
#line 282 "clargv.f"
			r__.r = z__1.r, r__.i = z__1.i;
#line 283 "clargv.f"
/* L40: */
#line 283 "clargv.f"
		    }
#line 284 "clargv.f"
		}
#line 285 "clargv.f"
	    }
#line 286 "clargv.f"
	}
#line 287 "clargv.f"
L50:
#line 288 "clargv.f"
	c__[ic] = cs;
#line 289 "clargv.f"
	i__2 = iy;
#line 289 "clargv.f"
	y[i__2].r = sn.r, y[i__2].i = sn.i;
#line 290 "clargv.f"
	i__2 = ix;
#line 290 "clargv.f"
	x[i__2].r = r__.r, x[i__2].i = r__.i;
#line 291 "clargv.f"
	ic += *incc;
#line 292 "clargv.f"
	iy += *incy;
#line 293 "clargv.f"
	ix += *incx;
#line 294 "clargv.f"
/* L60: */
#line 294 "clargv.f"
    }
#line 295 "clargv.f"
    return 0;

/*     End of CLARGV */

} /* clargv_ */

