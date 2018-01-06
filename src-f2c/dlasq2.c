#line 1 "dlasq2.f"
/* dlasq2.f -- translated by f2c (version 20100827).
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

#line 1 "dlasq2.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__10 = 10;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__11 = 11;

/* > \brief \b DLASQ2 computes all the eigenvalues of the symmetric positive definite tridiagonal matrix assoc
iated with the qd Array Z to high relative accuracy. Used by sbdsqr and sstegr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASQ2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASQ2( N, Z, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ2 computes all the eigenvalues of the symmetric positive */
/* > definite tridiagonal matrix associated with the qd array Z to high */
/* > relative accuracy are computed to high relative accuracy, in the */
/* > absence of denormalization, underflow and overflow. */
/* > */
/* > To see the relation of Z to the tridiagonal matrix, let L be a */
/* > unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and */
/* > let U be an upper bidiagonal matrix with 1's above and diagonal */
/* > Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the */
/* > symmetric tridiagonal to which it is similar. */
/* > */
/* > Note : DLASQ2 defines a logical variable, IEEE, which is true */
/* > on machines which follow ieee-754 floating-point standard in their */
/* > handling of infinities and NaNs, and false otherwise. This variable */
/* > is passed to DLASQ3. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >        The number of rows and columns in the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( 4*N ) */
/* >        On entry Z holds the qd array. On exit, entries 1 to N hold */
/* >        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the */
/* >        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If */
/* >        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 ) */
/* >        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of */
/* >        shifts that failed. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >        = 0: successful exit */
/* >        < 0: if the i-th argument is a scalar and had an illegal */
/* >             value, then INFO = -i, if the i-th argument is an */
/* >             array and the j-entry had an illegal value, then */
/* >             INFO = -(i*100+j) */
/* >        > 0: the algorithm failed */
/* >              = 1, a split was marked by a positive value in E */
/* >              = 2, current block of Z not diagonalized after 100*N */
/* >                   iterations (in inner while loop).  On exit Z holds */
/* >                   a qd array with the same eigenvalues as the given Z. */
/* >              = 3, termination criterion of outer while loop not met */
/* >                   (program created more than N unreduced blocks) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Local Variables: I0:N0 defines a current unreduced segment of Z. */
/* >  The shifts are accumulated in SIGMA. Iteration count is in ITER. */
/* >  Ping-pong is controlled by PP (alternates between 0 and 1). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlasq2_(integer *n, doublereal *z__, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__, e, g;
    static integer k;
    static doublereal s, t;
    static integer i0, i1, i4, n0, n1;
    static doublereal dn;
    static integer pp;
    static doublereal dn1, dn2, dee, eps, tau, tol;
    static integer ipn4;
    static doublereal tol2;
    static logical ieee;
    static integer nbig;
    static doublereal dmin__, emin, emax;
    static integer kmin, ndiv, iter;
    static doublereal qmin, temp, qmax, zmax;
    static integer splt;
    static doublereal dmin1, dmin2;
    static integer nfail;
    static doublereal desig, trace, sigma;
    static integer iinfo;
    static doublereal tempe, tempq;
    static integer ttype;
    extern /* Subroutine */ int dlasq3_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, logical *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal deemin;
    static integer iwhila, iwhilb;
    static doublereal oldemn, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlasrt_(char *, integer *, doublereal *, 
	    integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */
/*     (in case DLASQ2 is not called by DLASQ1) */

#line 162 "dlasq2.f"
    /* Parameter adjustments */
#line 162 "dlasq2.f"
    --z__;
#line 162 "dlasq2.f"

#line 162 "dlasq2.f"
    /* Function Body */
#line 162 "dlasq2.f"
    *info = 0;
#line 163 "dlasq2.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 164 "dlasq2.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 165 "dlasq2.f"
    tol = eps * 100.;
/* Computing 2nd power */
#line 166 "dlasq2.f"
    d__1 = tol;
#line 166 "dlasq2.f"
    tol2 = d__1 * d__1;

#line 168 "dlasq2.f"
    if (*n < 0) {
#line 169 "dlasq2.f"
	*info = -1;
#line 170 "dlasq2.f"
	xerbla_("DLASQ2", &c__1, (ftnlen)6);
#line 171 "dlasq2.f"
	return 0;
#line 172 "dlasq2.f"
    } else if (*n == 0) {
#line 173 "dlasq2.f"
	return 0;
#line 174 "dlasq2.f"
    } else if (*n == 1) {

/*        1-by-1 case. */

#line 178 "dlasq2.f"
	if (z__[1] < 0.) {
#line 179 "dlasq2.f"
	    *info = -201;
#line 180 "dlasq2.f"
	    xerbla_("DLASQ2", &c__2, (ftnlen)6);
#line 181 "dlasq2.f"
	}
#line 182 "dlasq2.f"
	return 0;
#line 183 "dlasq2.f"
    } else if (*n == 2) {

/*        2-by-2 case. */

#line 187 "dlasq2.f"
	if (z__[2] < 0. || z__[3] < 0.) {
#line 188 "dlasq2.f"
	    *info = -2;
#line 189 "dlasq2.f"
	    xerbla_("DLASQ2", &c__2, (ftnlen)6);
#line 190 "dlasq2.f"
	    return 0;
#line 191 "dlasq2.f"
	} else if (z__[3] > z__[1]) {
#line 192 "dlasq2.f"
	    d__ = z__[3];
#line 193 "dlasq2.f"
	    z__[3] = z__[1];
#line 194 "dlasq2.f"
	    z__[1] = d__;
#line 195 "dlasq2.f"
	}
#line 196 "dlasq2.f"
	z__[5] = z__[1] + z__[2] + z__[3];
#line 197 "dlasq2.f"
	if (z__[2] > z__[3] * tol2) {
#line 198 "dlasq2.f"
	    t = (z__[1] - z__[3] + z__[2]) * .5;
#line 199 "dlasq2.f"
	    s = z__[3] * (z__[2] / t);
#line 200 "dlasq2.f"
	    if (s <= t) {
#line 201 "dlasq2.f"
		s = z__[3] * (z__[2] / (t * (sqrt(s / t + 1.) + 1.)));
#line 202 "dlasq2.f"
	    } else {
#line 203 "dlasq2.f"
		s = z__[3] * (z__[2] / (t + sqrt(t) * sqrt(t + s)));
#line 204 "dlasq2.f"
	    }
#line 205 "dlasq2.f"
	    t = z__[1] + (s + z__[2]);
#line 206 "dlasq2.f"
	    z__[3] *= z__[1] / t;
#line 207 "dlasq2.f"
	    z__[1] = t;
#line 208 "dlasq2.f"
	}
#line 209 "dlasq2.f"
	z__[2] = z__[3];
#line 210 "dlasq2.f"
	z__[6] = z__[2] + z__[1];
#line 211 "dlasq2.f"
	return 0;
#line 212 "dlasq2.f"
    }

/*     Check for negative data and compute sums of q's and e's. */

#line 216 "dlasq2.f"
    z__[*n * 2] = 0.;
#line 217 "dlasq2.f"
    emin = z__[2];
#line 218 "dlasq2.f"
    qmax = 0.;
#line 219 "dlasq2.f"
    zmax = 0.;
#line 220 "dlasq2.f"
    d__ = 0.;
#line 221 "dlasq2.f"
    e = 0.;

#line 223 "dlasq2.f"
    i__1 = *n - 1 << 1;
#line 223 "dlasq2.f"
    for (k = 1; k <= i__1; k += 2) {
#line 224 "dlasq2.f"
	if (z__[k] < 0.) {
#line 225 "dlasq2.f"
	    *info = -(k + 200);
#line 226 "dlasq2.f"
	    xerbla_("DLASQ2", &c__2, (ftnlen)6);
#line 227 "dlasq2.f"
	    return 0;
#line 228 "dlasq2.f"
	} else if (z__[k + 1] < 0.) {
#line 229 "dlasq2.f"
	    *info = -(k + 201);
#line 230 "dlasq2.f"
	    xerbla_("DLASQ2", &c__2, (ftnlen)6);
#line 231 "dlasq2.f"
	    return 0;
#line 232 "dlasq2.f"
	}
#line 233 "dlasq2.f"
	d__ += z__[k];
#line 234 "dlasq2.f"
	e += z__[k + 1];
/* Computing MAX */
#line 235 "dlasq2.f"
	d__1 = qmax, d__2 = z__[k];
#line 235 "dlasq2.f"
	qmax = max(d__1,d__2);
/* Computing MIN */
#line 236 "dlasq2.f"
	d__1 = emin, d__2 = z__[k + 1];
#line 236 "dlasq2.f"
	emin = min(d__1,d__2);
/* Computing MAX */
#line 237 "dlasq2.f"
	d__1 = max(qmax,zmax), d__2 = z__[k + 1];
#line 237 "dlasq2.f"
	zmax = max(d__1,d__2);
#line 238 "dlasq2.f"
/* L10: */
#line 238 "dlasq2.f"
    }
#line 239 "dlasq2.f"
    if (z__[(*n << 1) - 1] < 0.) {
#line 240 "dlasq2.f"
	*info = -((*n << 1) + 199);
#line 241 "dlasq2.f"
	xerbla_("DLASQ2", &c__2, (ftnlen)6);
#line 242 "dlasq2.f"
	return 0;
#line 243 "dlasq2.f"
    }
#line 244 "dlasq2.f"
    d__ += z__[(*n << 1) - 1];
/* Computing MAX */
#line 245 "dlasq2.f"
    d__1 = qmax, d__2 = z__[(*n << 1) - 1];
#line 245 "dlasq2.f"
    qmax = max(d__1,d__2);
#line 246 "dlasq2.f"
    zmax = max(qmax,zmax);

/*     Check for diagonality. */

#line 250 "dlasq2.f"
    if (e == 0.) {
#line 251 "dlasq2.f"
	i__1 = *n;
#line 251 "dlasq2.f"
	for (k = 2; k <= i__1; ++k) {
#line 252 "dlasq2.f"
	    z__[k] = z__[(k << 1) - 1];
#line 253 "dlasq2.f"
/* L20: */
#line 253 "dlasq2.f"
	}
#line 254 "dlasq2.f"
	dlasrt_("D", n, &z__[1], &iinfo, (ftnlen)1);
#line 255 "dlasq2.f"
	z__[(*n << 1) - 1] = d__;
#line 256 "dlasq2.f"
	return 0;
#line 257 "dlasq2.f"
    }

#line 259 "dlasq2.f"
    trace = d__ + e;

/*     Check for zero data. */

#line 263 "dlasq2.f"
    if (trace == 0.) {
#line 264 "dlasq2.f"
	z__[(*n << 1) - 1] = 0.;
#line 265 "dlasq2.f"
	return 0;
#line 266 "dlasq2.f"
    }

/*     Check whether the machine is IEEE conformable. */

#line 270 "dlasq2.f"
    ieee = ilaenv_(&c__10, "DLASQ2", "N", &c__1, &c__2, &c__3, &c__4, (ftnlen)
	    6, (ftnlen)1) == 1 && ilaenv_(&c__11, "DLASQ2", "N", &c__1, &c__2,
	     &c__3, &c__4, (ftnlen)6, (ftnlen)1) == 1;

/*     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...). */

#line 275 "dlasq2.f"
    for (k = *n << 1; k >= 2; k += -2) {
#line 276 "dlasq2.f"
	z__[k * 2] = 0.;
#line 277 "dlasq2.f"
	z__[(k << 1) - 1] = z__[k];
#line 278 "dlasq2.f"
	z__[(k << 1) - 2] = 0.;
#line 279 "dlasq2.f"
	z__[(k << 1) - 3] = z__[k - 1];
#line 280 "dlasq2.f"
/* L30: */
#line 280 "dlasq2.f"
    }

#line 282 "dlasq2.f"
    i0 = 1;
#line 283 "dlasq2.f"
    n0 = *n;

/*     Reverse the qd-array, if warranted. */

#line 287 "dlasq2.f"
    if (z__[(i0 << 2) - 3] * 1.5 < z__[(n0 << 2) - 3]) {
#line 288 "dlasq2.f"
	ipn4 = i0 + n0 << 2;
#line 289 "dlasq2.f"
	i__1 = i0 + n0 - 1 << 1;
#line 289 "dlasq2.f"
	for (i4 = i0 << 2; i4 <= i__1; i4 += 4) {
#line 290 "dlasq2.f"
	    temp = z__[i4 - 3];
#line 291 "dlasq2.f"
	    z__[i4 - 3] = z__[ipn4 - i4 - 3];
#line 292 "dlasq2.f"
	    z__[ipn4 - i4 - 3] = temp;
#line 293 "dlasq2.f"
	    temp = z__[i4 - 1];
#line 294 "dlasq2.f"
	    z__[i4 - 1] = z__[ipn4 - i4 - 5];
#line 295 "dlasq2.f"
	    z__[ipn4 - i4 - 5] = temp;
#line 296 "dlasq2.f"
/* L40: */
#line 296 "dlasq2.f"
	}
#line 297 "dlasq2.f"
    }

/*     Initial split checking via dqd and Li's test. */

#line 301 "dlasq2.f"
    pp = 0;

#line 303 "dlasq2.f"
    for (k = 1; k <= 2; ++k) {

#line 305 "dlasq2.f"
	d__ = z__[(n0 << 2) + pp - 3];
#line 306 "dlasq2.f"
	i__1 = (i0 << 2) + pp;
#line 306 "dlasq2.f"
	for (i4 = (n0 - 1 << 2) + pp; i4 >= i__1; i4 += -4) {
#line 307 "dlasq2.f"
	    if (z__[i4 - 1] <= tol2 * d__) {
#line 308 "dlasq2.f"
		z__[i4 - 1] = -0.;
#line 309 "dlasq2.f"
		d__ = z__[i4 - 3];
#line 310 "dlasq2.f"
	    } else {
#line 311 "dlasq2.f"
		d__ = z__[i4 - 3] * (d__ / (d__ + z__[i4 - 1]));
#line 312 "dlasq2.f"
	    }
#line 313 "dlasq2.f"
/* L50: */
#line 313 "dlasq2.f"
	}

/*        dqd maps Z to ZZ plus Li's test. */

#line 317 "dlasq2.f"
	emin = z__[(i0 << 2) + pp + 1];
#line 318 "dlasq2.f"
	d__ = z__[(i0 << 2) + pp - 3];
#line 319 "dlasq2.f"
	i__1 = (n0 - 1 << 2) + pp;
#line 319 "dlasq2.f"
	for (i4 = (i0 << 2) + pp; i4 <= i__1; i4 += 4) {
#line 320 "dlasq2.f"
	    z__[i4 - (pp << 1) - 2] = d__ + z__[i4 - 1];
#line 321 "dlasq2.f"
	    if (z__[i4 - 1] <= tol2 * d__) {
#line 322 "dlasq2.f"
		z__[i4 - 1] = -0.;
#line 323 "dlasq2.f"
		z__[i4 - (pp << 1) - 2] = d__;
#line 324 "dlasq2.f"
		z__[i4 - (pp << 1)] = 0.;
#line 325 "dlasq2.f"
		d__ = z__[i4 + 1];
#line 326 "dlasq2.f"
	    } else if (safmin * z__[i4 + 1] < z__[i4 - (pp << 1) - 2] && 
		    safmin * z__[i4 - (pp << 1) - 2] < z__[i4 + 1]) {
#line 328 "dlasq2.f"
		temp = z__[i4 + 1] / z__[i4 - (pp << 1) - 2];
#line 329 "dlasq2.f"
		z__[i4 - (pp << 1)] = z__[i4 - 1] * temp;
#line 330 "dlasq2.f"
		d__ *= temp;
#line 331 "dlasq2.f"
	    } else {
#line 332 "dlasq2.f"
		z__[i4 - (pp << 1)] = z__[i4 + 1] * (z__[i4 - 1] / z__[i4 - (
			pp << 1) - 2]);
#line 333 "dlasq2.f"
		d__ = z__[i4 + 1] * (d__ / z__[i4 - (pp << 1) - 2]);
#line 334 "dlasq2.f"
	    }
/* Computing MIN */
#line 335 "dlasq2.f"
	    d__1 = emin, d__2 = z__[i4 - (pp << 1)];
#line 335 "dlasq2.f"
	    emin = min(d__1,d__2);
#line 336 "dlasq2.f"
/* L60: */
#line 336 "dlasq2.f"
	}
#line 337 "dlasq2.f"
	z__[(n0 << 2) - pp - 2] = d__;

/*        Now find qmax. */

#line 341 "dlasq2.f"
	qmax = z__[(i0 << 2) - pp - 2];
#line 342 "dlasq2.f"
	i__1 = (n0 << 2) - pp - 2;
#line 342 "dlasq2.f"
	for (i4 = (i0 << 2) - pp + 2; i4 <= i__1; i4 += 4) {
/* Computing MAX */
#line 343 "dlasq2.f"
	    d__1 = qmax, d__2 = z__[i4];
#line 343 "dlasq2.f"
	    qmax = max(d__1,d__2);
#line 344 "dlasq2.f"
/* L70: */
#line 344 "dlasq2.f"
	}

/*        Prepare for the next iteration on K. */

#line 348 "dlasq2.f"
	pp = 1 - pp;
#line 349 "dlasq2.f"
/* L80: */
#line 349 "dlasq2.f"
    }

/*     Initialise variables to pass to DLASQ3. */

#line 353 "dlasq2.f"
    ttype = 0;
#line 354 "dlasq2.f"
    dmin1 = 0.;
#line 355 "dlasq2.f"
    dmin2 = 0.;
#line 356 "dlasq2.f"
    dn = 0.;
#line 357 "dlasq2.f"
    dn1 = 0.;
#line 358 "dlasq2.f"
    dn2 = 0.;
#line 359 "dlasq2.f"
    g = 0.;
#line 360 "dlasq2.f"
    tau = 0.;

#line 362 "dlasq2.f"
    iter = 2;
#line 363 "dlasq2.f"
    nfail = 0;
#line 364 "dlasq2.f"
    ndiv = n0 - i0 << 1;

#line 366 "dlasq2.f"
    i__1 = *n + 1;
#line 366 "dlasq2.f"
    for (iwhila = 1; iwhila <= i__1; ++iwhila) {
#line 367 "dlasq2.f"
	if (n0 < 1) {
#line 367 "dlasq2.f"
	    goto L170;
#line 367 "dlasq2.f"
	}

/*        While array unfinished do */

/*        E(N0) holds the value of SIGMA when submatrix in I0:N0 */
/*        splits from the rest of the array, but is negated. */

#line 375 "dlasq2.f"
	desig = 0.;
#line 376 "dlasq2.f"
	if (n0 == *n) {
#line 377 "dlasq2.f"
	    sigma = 0.;
#line 378 "dlasq2.f"
	} else {
#line 379 "dlasq2.f"
	    sigma = -z__[(n0 << 2) - 1];
#line 380 "dlasq2.f"
	}
#line 381 "dlasq2.f"
	if (sigma < 0.) {
#line 382 "dlasq2.f"
	    *info = 1;
#line 383 "dlasq2.f"
	    return 0;
#line 384 "dlasq2.f"
	}

/*        Find last unreduced submatrix's top index I0, find QMAX and */
/*        EMIN. Find Gershgorin-type bound if Q's much greater than E's. */

#line 389 "dlasq2.f"
	emax = 0.;
#line 390 "dlasq2.f"
	if (n0 > i0) {
#line 391 "dlasq2.f"
	    emin = (d__1 = z__[(n0 << 2) - 5], abs(d__1));
#line 392 "dlasq2.f"
	} else {
#line 393 "dlasq2.f"
	    emin = 0.;
#line 394 "dlasq2.f"
	}
#line 395 "dlasq2.f"
	qmin = z__[(n0 << 2) - 3];
#line 396 "dlasq2.f"
	qmax = qmin;
#line 397 "dlasq2.f"
	for (i4 = n0 << 2; i4 >= 8; i4 += -4) {
#line 398 "dlasq2.f"
	    if (z__[i4 - 5] <= 0.) {
#line 398 "dlasq2.f"
		goto L100;
#line 398 "dlasq2.f"
	    }
#line 400 "dlasq2.f"
	    if (qmin >= emax * 4.) {
/* Computing MIN */
#line 401 "dlasq2.f"
		d__1 = qmin, d__2 = z__[i4 - 3];
#line 401 "dlasq2.f"
		qmin = min(d__1,d__2);
/* Computing MAX */
#line 402 "dlasq2.f"
		d__1 = emax, d__2 = z__[i4 - 5];
#line 402 "dlasq2.f"
		emax = max(d__1,d__2);
#line 403 "dlasq2.f"
	    }
/* Computing MAX */
#line 404 "dlasq2.f"
	    d__1 = qmax, d__2 = z__[i4 - 7] + z__[i4 - 5];
#line 404 "dlasq2.f"
	    qmax = max(d__1,d__2);
/* Computing MIN */
#line 405 "dlasq2.f"
	    d__1 = emin, d__2 = z__[i4 - 5];
#line 405 "dlasq2.f"
	    emin = min(d__1,d__2);
#line 406 "dlasq2.f"
/* L90: */
#line 406 "dlasq2.f"
	}
#line 407 "dlasq2.f"
	i4 = 4;

#line 409 "dlasq2.f"
L100:
#line 410 "dlasq2.f"
	i0 = i4 / 4;
#line 411 "dlasq2.f"
	pp = 0;

#line 413 "dlasq2.f"
	if (n0 - i0 > 1) {
#line 414 "dlasq2.f"
	    dee = z__[(i0 << 2) - 3];
#line 415 "dlasq2.f"
	    deemin = dee;
#line 416 "dlasq2.f"
	    kmin = i0;
#line 417 "dlasq2.f"
	    i__2 = (n0 << 2) - 3;
#line 417 "dlasq2.f"
	    for (i4 = (i0 << 2) + 1; i4 <= i__2; i4 += 4) {
#line 418 "dlasq2.f"
		dee = z__[i4] * (dee / (dee + z__[i4 - 2]));
#line 419 "dlasq2.f"
		if (dee <= deemin) {
#line 420 "dlasq2.f"
		    deemin = dee;
#line 421 "dlasq2.f"
		    kmin = (i4 + 3) / 4;
#line 422 "dlasq2.f"
		}
#line 423 "dlasq2.f"
/* L110: */
#line 423 "dlasq2.f"
	    }
#line 424 "dlasq2.f"
	    if (kmin - i0 << 1 < n0 - kmin && deemin <= z__[(n0 << 2) - 3] * 
		    .5) {
#line 426 "dlasq2.f"
		ipn4 = i0 + n0 << 2;
#line 427 "dlasq2.f"
		pp = 2;
#line 428 "dlasq2.f"
		i__2 = i0 + n0 - 1 << 1;
#line 428 "dlasq2.f"
		for (i4 = i0 << 2; i4 <= i__2; i4 += 4) {
#line 429 "dlasq2.f"
		    temp = z__[i4 - 3];
#line 430 "dlasq2.f"
		    z__[i4 - 3] = z__[ipn4 - i4 - 3];
#line 431 "dlasq2.f"
		    z__[ipn4 - i4 - 3] = temp;
#line 432 "dlasq2.f"
		    temp = z__[i4 - 2];
#line 433 "dlasq2.f"
		    z__[i4 - 2] = z__[ipn4 - i4 - 2];
#line 434 "dlasq2.f"
		    z__[ipn4 - i4 - 2] = temp;
#line 435 "dlasq2.f"
		    temp = z__[i4 - 1];
#line 436 "dlasq2.f"
		    z__[i4 - 1] = z__[ipn4 - i4 - 5];
#line 437 "dlasq2.f"
		    z__[ipn4 - i4 - 5] = temp;
#line 438 "dlasq2.f"
		    temp = z__[i4];
#line 439 "dlasq2.f"
		    z__[i4] = z__[ipn4 - i4 - 4];
#line 440 "dlasq2.f"
		    z__[ipn4 - i4 - 4] = temp;
#line 441 "dlasq2.f"
/* L120: */
#line 441 "dlasq2.f"
		}
#line 442 "dlasq2.f"
	    }
#line 443 "dlasq2.f"
	}

/*        Put -(initial shift) into DMIN. */

/* Computing MAX */
#line 447 "dlasq2.f"
	d__1 = 0., d__2 = qmin - sqrt(qmin) * 2. * sqrt(emax);
#line 447 "dlasq2.f"
	dmin__ = -max(d__1,d__2);

/*        Now I0:N0 is unreduced. */
/*        PP = 0 for ping, PP = 1 for pong. */
/*        PP = 2 indicates that flipping was applied to the Z array and */
/*               and that the tests for deflation upon entry in DLASQ3 */
/*               should not be performed. */

#line 455 "dlasq2.f"
	nbig = (n0 - i0 + 1) * 100;
#line 456 "dlasq2.f"
	i__2 = nbig;
#line 456 "dlasq2.f"
	for (iwhilb = 1; iwhilb <= i__2; ++iwhilb) {
#line 457 "dlasq2.f"
	    if (i0 > n0) {
#line 457 "dlasq2.f"
		goto L150;
#line 457 "dlasq2.f"
	    }

/*           While submatrix unfinished take a good dqds step. */

#line 462 "dlasq2.f"
	    dlasq3_(&i0, &n0, &z__[1], &pp, &dmin__, &sigma, &desig, &qmax, &
		    nfail, &iter, &ndiv, &ieee, &ttype, &dmin1, &dmin2, &dn, &
		    dn1, &dn2, &g, &tau);

#line 466 "dlasq2.f"
	    pp = 1 - pp;

/*           When EMIN is very small check for splits. */

#line 470 "dlasq2.f"
	    if (pp == 0 && n0 - i0 >= 3) {
#line 471 "dlasq2.f"
		if (z__[n0 * 4] <= tol2 * qmax || z__[(n0 << 2) - 1] <= tol2 *
			 sigma) {
#line 473 "dlasq2.f"
		    splt = i0 - 1;
#line 474 "dlasq2.f"
		    qmax = z__[(i0 << 2) - 3];
#line 475 "dlasq2.f"
		    emin = z__[(i0 << 2) - 1];
#line 476 "dlasq2.f"
		    oldemn = z__[i0 * 4];
#line 477 "dlasq2.f"
		    i__3 = n0 - 3 << 2;
#line 477 "dlasq2.f"
		    for (i4 = i0 << 2; i4 <= i__3; i4 += 4) {
#line 478 "dlasq2.f"
			if (z__[i4] <= tol2 * z__[i4 - 3] || z__[i4 - 1] <= 
				tol2 * sigma) {
#line 480 "dlasq2.f"
			    z__[i4 - 1] = -sigma;
#line 481 "dlasq2.f"
			    splt = i4 / 4;
#line 482 "dlasq2.f"
			    qmax = 0.;
#line 483 "dlasq2.f"
			    emin = z__[i4 + 3];
#line 484 "dlasq2.f"
			    oldemn = z__[i4 + 4];
#line 485 "dlasq2.f"
			} else {
/* Computing MAX */
#line 486 "dlasq2.f"
			    d__1 = qmax, d__2 = z__[i4 + 1];
#line 486 "dlasq2.f"
			    qmax = max(d__1,d__2);
/* Computing MIN */
#line 487 "dlasq2.f"
			    d__1 = emin, d__2 = z__[i4 - 1];
#line 487 "dlasq2.f"
			    emin = min(d__1,d__2);
/* Computing MIN */
#line 488 "dlasq2.f"
			    d__1 = oldemn, d__2 = z__[i4];
#line 488 "dlasq2.f"
			    oldemn = min(d__1,d__2);
#line 489 "dlasq2.f"
			}
#line 490 "dlasq2.f"
/* L130: */
#line 490 "dlasq2.f"
		    }
#line 491 "dlasq2.f"
		    z__[(n0 << 2) - 1] = emin;
#line 492 "dlasq2.f"
		    z__[n0 * 4] = oldemn;
#line 493 "dlasq2.f"
		    i0 = splt + 1;
#line 494 "dlasq2.f"
		}
#line 495 "dlasq2.f"
	    }

#line 497 "dlasq2.f"
/* L140: */
#line 497 "dlasq2.f"
	}

#line 499 "dlasq2.f"
	*info = 2;

/*        Maximum number of iterations exceeded, restore the shift */
/*        SIGMA and place the new d's and e's in a qd array. */
/*        This might need to be done for several blocks */

#line 505 "dlasq2.f"
	i1 = i0;
#line 506 "dlasq2.f"
	n1 = n0;
#line 507 "dlasq2.f"
L145:
#line 508 "dlasq2.f"
	tempq = z__[(i0 << 2) - 3];
#line 509 "dlasq2.f"
	z__[(i0 << 2) - 3] += sigma;
#line 510 "dlasq2.f"
	i__2 = n0;
#line 510 "dlasq2.f"
	for (k = i0 + 1; k <= i__2; ++k) {
#line 511 "dlasq2.f"
	    tempe = z__[(k << 2) - 5];
#line 512 "dlasq2.f"
	    z__[(k << 2) - 5] *= tempq / z__[(k << 2) - 7];
#line 513 "dlasq2.f"
	    tempq = z__[(k << 2) - 3];
#line 514 "dlasq2.f"
	    z__[(k << 2) - 3] = z__[(k << 2) - 3] + sigma + tempe - z__[(k << 
		    2) - 5];
#line 515 "dlasq2.f"
	}

/*        Prepare to do this on the previous block if there is one */

#line 519 "dlasq2.f"
	if (i1 > 1) {
#line 520 "dlasq2.f"
	    n1 = i1 - 1;
#line 521 "dlasq2.f"
	    while(i1 >= 2 && z__[(i1 << 2) - 5] >= 0.) {
#line 522 "dlasq2.f"
		--i1;
#line 523 "dlasq2.f"
	    }
#line 524 "dlasq2.f"
	    sigma = -z__[(n1 << 2) - 1];
#line 525 "dlasq2.f"
	    goto L145;
#line 526 "dlasq2.f"
	}
#line 528 "dlasq2.f"
	i__2 = *n;
#line 528 "dlasq2.f"
	for (k = 1; k <= i__2; ++k) {
#line 529 "dlasq2.f"
	    z__[(k << 1) - 1] = z__[(k << 2) - 3];

/*        Only the block 1..N0 is unfinished.  The rest of the e's */
/*        must be essentially zero, although sometimes other data */
/*        has been stored in them. */

#line 535 "dlasq2.f"
	    if (k < n0) {
#line 536 "dlasq2.f"
		z__[k * 2] = z__[(k << 2) - 1];
#line 537 "dlasq2.f"
	    } else {
#line 538 "dlasq2.f"
		z__[k * 2] = 0.;
#line 539 "dlasq2.f"
	    }
#line 540 "dlasq2.f"
	}
#line 541 "dlasq2.f"
	return 0;

/*        end IWHILB */

#line 545 "dlasq2.f"
L150:

#line 547 "dlasq2.f"
/* L160: */
#line 547 "dlasq2.f"
	;
#line 547 "dlasq2.f"
    }

#line 549 "dlasq2.f"
    *info = 3;
#line 550 "dlasq2.f"
    return 0;

/*     end IWHILA */

#line 554 "dlasq2.f"
L170:

/*     Move q's to the front. */

#line 558 "dlasq2.f"
    i__1 = *n;
#line 558 "dlasq2.f"
    for (k = 2; k <= i__1; ++k) {
#line 559 "dlasq2.f"
	z__[k] = z__[(k << 2) - 3];
#line 560 "dlasq2.f"
/* L180: */
#line 560 "dlasq2.f"
    }

/*     Sort and compute sum of eigenvalues. */

#line 564 "dlasq2.f"
    dlasrt_("D", n, &z__[1], &iinfo, (ftnlen)1);

#line 566 "dlasq2.f"
    e = 0.;
#line 567 "dlasq2.f"
    for (k = *n; k >= 1; --k) {
#line 568 "dlasq2.f"
	e += z__[k];
#line 569 "dlasq2.f"
/* L190: */
#line 569 "dlasq2.f"
    }

/*     Store trace, sum(eigenvalues) and information on performance. */

#line 573 "dlasq2.f"
    z__[(*n << 1) + 1] = trace;
#line 574 "dlasq2.f"
    z__[(*n << 1) + 2] = e;
#line 575 "dlasq2.f"
    z__[(*n << 1) + 3] = (doublereal) iter;
/* Computing 2nd power */
#line 576 "dlasq2.f"
    i__1 = *n;
#line 576 "dlasq2.f"
    z__[(*n << 1) + 4] = (doublereal) ndiv / (doublereal) (i__1 * i__1);
#line 577 "dlasq2.f"
    z__[(*n << 1) + 5] = nfail * 100. / (doublereal) iter;
#line 578 "dlasq2.f"
    return 0;

/*     End of DLASQ2 */

} /* dlasq2_ */

