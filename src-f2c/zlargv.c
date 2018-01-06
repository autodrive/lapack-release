#line 1 "zlargv.f"
/* zlargv.f -- translated by f2c (version 20100827).
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

#line 1 "zlargv.f"
/* > \brief \b ZLARGV generates a vector of plane rotations with real cosines and complex sines. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLARGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlargv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlargv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlargv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARGV( N, X, INCX, Y, INCY, C, INCC ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCC, INCX, INCY, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( * ) */
/*       COMPLEX*16         X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARGV generates a vector of complex plane rotations with real */
/* > cosines, determined by elements of the complex vectors x and y. */
/* > For i = 1,2,...,n */
/* > */
/* >    (        c(i)   s(i) ) ( x(i) ) = ( r(i) ) */
/* >    ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  ) */
/* > */
/* >    where c(i)**2 + ABS(s(i))**2 = 1 */
/* > */
/* > The following conventions are used (these are the same as in ZLARTG, */
/* > but differ from the BLAS1 routine ZROTG): */
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
/* >          X is COMPLEX*16 array, dimension (1+(N-1)*INCX) */
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
/* >          Y is COMPLEX*16 array, dimension (1+(N-1)*INCY) */
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
/* >          C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

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
/* Subroutine */ int zlargv_(integer *n, doublecomplex *x, integer *incx, 
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
    static doublereal safmn2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal safmx2;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;


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
#line 179 "zlargv.f"
    /* Parameter adjustments */
#line 179 "zlargv.f"
    --c__;
#line 179 "zlargv.f"
    --y;
#line 179 "zlargv.f"
    --x;
#line 179 "zlargv.f"

#line 179 "zlargv.f"
    /* Function Body */
#line 179 "zlargv.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 180 "zlargv.f"
    eps = dlamch_("E", (ftnlen)1);
#line 181 "zlargv.f"
    d__1 = dlamch_("B", (ftnlen)1);
#line 181 "zlargv.f"
    i__1 = (integer) (log(safmin / eps) / log(dlamch_("B", (ftnlen)1)) / 2.);
#line 181 "zlargv.f"
    safmn2 = pow_di(&d__1, &i__1);
#line 183 "zlargv.f"
    safmx2 = 1. / safmn2;
/*     END IF */
#line 185 "zlargv.f"
    ix = 1;
#line 186 "zlargv.f"
    iy = 1;
#line 187 "zlargv.f"
    ic = 1;
#line 188 "zlargv.f"
    i__1 = *n;
#line 188 "zlargv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 189 "zlargv.f"
	i__2 = ix;
#line 189 "zlargv.f"
	f.r = x[i__2].r, f.i = x[i__2].i;
#line 190 "zlargv.f"
	i__2 = iy;
#line 190 "zlargv.f"
	g.r = y[i__2].r, g.i = y[i__2].i;

/*        Use identical algorithm as in ZLARTG */

/* Computing MAX */
/* Computing MAX */
#line 194 "zlargv.f"
	d__7 = (d__1 = f.r, abs(d__1)), d__8 = (d__2 = d_imag(&f), abs(d__2));
/* Computing MAX */
#line 194 "zlargv.f"
	d__9 = (d__3 = g.r, abs(d__3)), d__10 = (d__4 = d_imag(&g), abs(d__4))
		;
#line 194 "zlargv.f"
	d__5 = max(d__7,d__8), d__6 = max(d__9,d__10);
#line 194 "zlargv.f"
	scale = max(d__5,d__6);
#line 195 "zlargv.f"
	fs.r = f.r, fs.i = f.i;
#line 196 "zlargv.f"
	gs.r = g.r, gs.i = g.i;
#line 197 "zlargv.f"
	count = 0;
#line 198 "zlargv.f"
	if (scale >= safmx2) {
#line 199 "zlargv.f"
L10:
#line 200 "zlargv.f"
	    ++count;
#line 201 "zlargv.f"
	    z__1.r = safmn2 * fs.r, z__1.i = safmn2 * fs.i;
#line 201 "zlargv.f"
	    fs.r = z__1.r, fs.i = z__1.i;
#line 202 "zlargv.f"
	    z__1.r = safmn2 * gs.r, z__1.i = safmn2 * gs.i;
#line 202 "zlargv.f"
	    gs.r = z__1.r, gs.i = z__1.i;
#line 203 "zlargv.f"
	    scale *= safmn2;
#line 204 "zlargv.f"
	    if (scale >= safmx2) {
#line 204 "zlargv.f"
		goto L10;
#line 204 "zlargv.f"
	    }
#line 206 "zlargv.f"
	} else if (scale <= safmn2) {
#line 207 "zlargv.f"
	    if (g.r == 0. && g.i == 0.) {
#line 208 "zlargv.f"
		cs = 1.;
#line 209 "zlargv.f"
		sn.r = 0., sn.i = 0.;
#line 210 "zlargv.f"
		r__.r = f.r, r__.i = f.i;
#line 211 "zlargv.f"
		goto L50;
#line 212 "zlargv.f"
	    }
#line 213 "zlargv.f"
L20:
#line 214 "zlargv.f"
	    --count;
#line 215 "zlargv.f"
	    z__1.r = safmx2 * fs.r, z__1.i = safmx2 * fs.i;
#line 215 "zlargv.f"
	    fs.r = z__1.r, fs.i = z__1.i;
#line 216 "zlargv.f"
	    z__1.r = safmx2 * gs.r, z__1.i = safmx2 * gs.i;
#line 216 "zlargv.f"
	    gs.r = z__1.r, gs.i = z__1.i;
#line 217 "zlargv.f"
	    scale *= safmx2;
#line 218 "zlargv.f"
	    if (scale <= safmn2) {
#line 218 "zlargv.f"
		goto L20;
#line 218 "zlargv.f"
	    }
#line 220 "zlargv.f"
	}
/* Computing 2nd power */
#line 221 "zlargv.f"
	d__1 = fs.r;
/* Computing 2nd power */
#line 221 "zlargv.f"
	d__2 = d_imag(&fs);
#line 221 "zlargv.f"
	f2 = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
#line 222 "zlargv.f"
	d__1 = gs.r;
/* Computing 2nd power */
#line 222 "zlargv.f"
	d__2 = d_imag(&gs);
#line 222 "zlargv.f"
	g2 = d__1 * d__1 + d__2 * d__2;
#line 223 "zlargv.f"
	if (f2 <= max(g2,1.) * safmin) {

/*           This is a rare case: F is very small. */

#line 227 "zlargv.f"
	    if (f.r == 0. && f.i == 0.) {
#line 228 "zlargv.f"
		cs = 0.;
#line 229 "zlargv.f"
		d__2 = g.r;
#line 229 "zlargv.f"
		d__3 = d_imag(&g);
#line 229 "zlargv.f"
		d__1 = dlapy2_(&d__2, &d__3);
#line 229 "zlargv.f"
		r__.r = d__1, r__.i = 0.;
/*              Do complex/real division explicitly with two real */
/*              divisions */
#line 232 "zlargv.f"
		d__1 = gs.r;
#line 232 "zlargv.f"
		d__2 = d_imag(&gs);
#line 232 "zlargv.f"
		d__ = dlapy2_(&d__1, &d__2);
#line 233 "zlargv.f"
		d__1 = gs.r / d__;
#line 233 "zlargv.f"
		d__2 = -d_imag(&gs) / d__;
#line 233 "zlargv.f"
		z__1.r = d__1, z__1.i = d__2;
#line 233 "zlargv.f"
		sn.r = z__1.r, sn.i = z__1.i;
#line 234 "zlargv.f"
		goto L50;
#line 235 "zlargv.f"
	    }
#line 236 "zlargv.f"
	    d__1 = fs.r;
#line 236 "zlargv.f"
	    d__2 = d_imag(&fs);
#line 236 "zlargv.f"
	    f2s = dlapy2_(&d__1, &d__2);
/*           G2 and G2S are accurate */
/*           G2 is at least SAFMIN, and G2S is at least SAFMN2 */
#line 239 "zlargv.f"
	    g2s = sqrt(g2);
/*           Error in CS from underflow in F2S is at most */
/*           UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS */
/*           If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN, */
/*           and so CS .lt. sqrt(SAFMIN) */
/*           If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN */
/*           and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS) */
/*           Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S */
#line 247 "zlargv.f"
	    cs = f2s / g2s;
/*           Make sure abs(FF) = 1 */
/*           Do complex/real division explicitly with 2 real divisions */
/* Computing MAX */
#line 250 "zlargv.f"
	    d__3 = (d__1 = f.r, abs(d__1)), d__4 = (d__2 = d_imag(&f), abs(
		    d__2));
#line 250 "zlargv.f"
	    if (max(d__3,d__4) > 1.) {
#line 251 "zlargv.f"
		d__1 = f.r;
#line 251 "zlargv.f"
		d__2 = d_imag(&f);
#line 251 "zlargv.f"
		d__ = dlapy2_(&d__1, &d__2);
#line 252 "zlargv.f"
		d__1 = f.r / d__;
#line 252 "zlargv.f"
		d__2 = d_imag(&f) / d__;
#line 252 "zlargv.f"
		z__1.r = d__1, z__1.i = d__2;
#line 252 "zlargv.f"
		ff.r = z__1.r, ff.i = z__1.i;
#line 253 "zlargv.f"
	    } else {
#line 254 "zlargv.f"
		dr = safmx2 * f.r;
#line 255 "zlargv.f"
		di = safmx2 * d_imag(&f);
#line 256 "zlargv.f"
		d__ = dlapy2_(&dr, &di);
#line 257 "zlargv.f"
		d__1 = dr / d__;
#line 257 "zlargv.f"
		d__2 = di / d__;
#line 257 "zlargv.f"
		z__1.r = d__1, z__1.i = d__2;
#line 257 "zlargv.f"
		ff.r = z__1.r, ff.i = z__1.i;
#line 258 "zlargv.f"
	    }
#line 259 "zlargv.f"
	    d__1 = gs.r / g2s;
#line 259 "zlargv.f"
	    d__2 = -d_imag(&gs) / g2s;
#line 259 "zlargv.f"
	    z__2.r = d__1, z__2.i = d__2;
#line 259 "zlargv.f"
	    z__1.r = ff.r * z__2.r - ff.i * z__2.i, z__1.i = ff.r * z__2.i + 
		    ff.i * z__2.r;
#line 259 "zlargv.f"
	    sn.r = z__1.r, sn.i = z__1.i;
#line 260 "zlargv.f"
	    z__2.r = cs * f.r, z__2.i = cs * f.i;
#line 260 "zlargv.f"
	    z__3.r = sn.r * g.r - sn.i * g.i, z__3.i = sn.r * g.i + sn.i * 
		    g.r;
#line 260 "zlargv.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 260 "zlargv.f"
	    r__.r = z__1.r, r__.i = z__1.i;
#line 261 "zlargv.f"
	} else {

/*           This is the most common case. */
/*           Neither F2 nor F2/G2 are less than SAFMIN */
/*           F2S cannot overflow, and it is accurate */

#line 267 "zlargv.f"
	    f2s = sqrt(g2 / f2 + 1.);
/*           Do the F2S(real)*FS(complex) multiply with two real */
/*           multiplies */
#line 270 "zlargv.f"
	    d__1 = f2s * fs.r;
#line 270 "zlargv.f"
	    d__2 = f2s * d_imag(&fs);
#line 270 "zlargv.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 270 "zlargv.f"
	    r__.r = z__1.r, r__.i = z__1.i;
#line 271 "zlargv.f"
	    cs = 1. / f2s;
#line 272 "zlargv.f"
	    d__ = f2 + g2;
/*           Do complex/real division explicitly with two real divisions */
#line 274 "zlargv.f"
	    d__1 = r__.r / d__;
#line 274 "zlargv.f"
	    d__2 = d_imag(&r__) / d__;
#line 274 "zlargv.f"
	    z__1.r = d__1, z__1.i = d__2;
#line 274 "zlargv.f"
	    sn.r = z__1.r, sn.i = z__1.i;
#line 275 "zlargv.f"
	    d_cnjg(&z__2, &gs);
#line 275 "zlargv.f"
	    z__1.r = sn.r * z__2.r - sn.i * z__2.i, z__1.i = sn.r * z__2.i + 
		    sn.i * z__2.r;
#line 275 "zlargv.f"
	    sn.r = z__1.r, sn.i = z__1.i;
#line 276 "zlargv.f"
	    if (count != 0) {
#line 277 "zlargv.f"
		if (count > 0) {
#line 278 "zlargv.f"
		    i__2 = count;
#line 278 "zlargv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 279 "zlargv.f"
			z__1.r = safmx2 * r__.r, z__1.i = safmx2 * r__.i;
#line 279 "zlargv.f"
			r__.r = z__1.r, r__.i = z__1.i;
#line 280 "zlargv.f"
/* L30: */
#line 280 "zlargv.f"
		    }
#line 281 "zlargv.f"
		} else {
#line 282 "zlargv.f"
		    i__2 = -count;
#line 282 "zlargv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 283 "zlargv.f"
			z__1.r = safmn2 * r__.r, z__1.i = safmn2 * r__.i;
#line 283 "zlargv.f"
			r__.r = z__1.r, r__.i = z__1.i;
#line 284 "zlargv.f"
/* L40: */
#line 284 "zlargv.f"
		    }
#line 285 "zlargv.f"
		}
#line 286 "zlargv.f"
	    }
#line 287 "zlargv.f"
	}
#line 288 "zlargv.f"
L50:
#line 289 "zlargv.f"
	c__[ic] = cs;
#line 290 "zlargv.f"
	i__2 = iy;
#line 290 "zlargv.f"
	y[i__2].r = sn.r, y[i__2].i = sn.i;
#line 291 "zlargv.f"
	i__2 = ix;
#line 291 "zlargv.f"
	x[i__2].r = r__.r, x[i__2].i = r__.i;
#line 292 "zlargv.f"
	ic += *incc;
#line 293 "zlargv.f"
	iy += *incy;
#line 294 "zlargv.f"
	ix += *incx;
#line 295 "zlargv.f"
/* L60: */
#line 295 "zlargv.f"
    }
#line 296 "zlargv.f"
    return 0;

/*     End of ZLARGV */

} /* zlargv_ */

