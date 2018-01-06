#line 1 "slasq2.f"
/* slasq2.f -- translated by f2c (version 20100827).
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

#line 1 "slasq2.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b SLASQ2 computes all the eigenvalues of the symmetric positive definite tridiagonal matrix assoc
iated with the qd Array Z to high relative accuracy. Used by sbdsqr and sstegr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASQ2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASQ2( N, Z, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASQ2 computes all the eigenvalues of the symmetric positive */
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
/* > Note : SLASQ2 defines a logical variable, IEEE, which is true */
/* > on machines which follow ieee-754 floating-point standard in their */
/* > handling of infinities and NaNs, and false otherwise. This variable */
/* > is passed to SLASQ3. */
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
/* >          Z is REAL array, dimension ( 4*N ) */
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
/* Subroutine */ int slasq2_(integer *n, doublereal *z__, integer *info)
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
    extern /* Subroutine */ int slasq3_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, logical *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal deemin;
    extern doublereal slamch_(char *, ftnlen);
    static integer iwhila, iwhilb;
    static doublereal oldemn, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slasrt_(
	    char *, integer *, doublereal *, integer *, ftnlen);


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
/*     (in case SLASQ2 is not called by SLASQ1) */

#line 161 "slasq2.f"
    /* Parameter adjustments */
#line 161 "slasq2.f"
    --z__;
#line 161 "slasq2.f"

#line 161 "slasq2.f"
    /* Function Body */
#line 161 "slasq2.f"
    *info = 0;
#line 162 "slasq2.f"
    eps = slamch_("Precision", (ftnlen)9);
#line 163 "slasq2.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 164 "slasq2.f"
    tol = eps * 100.;
/* Computing 2nd power */
#line 165 "slasq2.f"
    d__1 = tol;
#line 165 "slasq2.f"
    tol2 = d__1 * d__1;

#line 167 "slasq2.f"
    if (*n < 0) {
#line 168 "slasq2.f"
	*info = -1;
#line 169 "slasq2.f"
	xerbla_("SLASQ2", &c__1, (ftnlen)6);
#line 170 "slasq2.f"
	return 0;
#line 171 "slasq2.f"
    } else if (*n == 0) {
#line 172 "slasq2.f"
	return 0;
#line 173 "slasq2.f"
    } else if (*n == 1) {

/*        1-by-1 case. */

#line 177 "slasq2.f"
	if (z__[1] < 0.) {
#line 178 "slasq2.f"
	    *info = -201;
#line 179 "slasq2.f"
	    xerbla_("SLASQ2", &c__2, (ftnlen)6);
#line 180 "slasq2.f"
	}
#line 181 "slasq2.f"
	return 0;
#line 182 "slasq2.f"
    } else if (*n == 2) {

/*        2-by-2 case. */

#line 186 "slasq2.f"
	if (z__[2] < 0. || z__[3] < 0.) {
#line 187 "slasq2.f"
	    *info = -2;
#line 188 "slasq2.f"
	    xerbla_("SLASQ2", &c__2, (ftnlen)6);
#line 189 "slasq2.f"
	    return 0;
#line 190 "slasq2.f"
	} else if (z__[3] > z__[1]) {
#line 191 "slasq2.f"
	    d__ = z__[3];
#line 192 "slasq2.f"
	    z__[3] = z__[1];
#line 193 "slasq2.f"
	    z__[1] = d__;
#line 194 "slasq2.f"
	}
#line 195 "slasq2.f"
	z__[5] = z__[1] + z__[2] + z__[3];
#line 196 "slasq2.f"
	if (z__[2] > z__[3] * tol2) {
#line 197 "slasq2.f"
	    t = (z__[1] - z__[3] + z__[2]) * .5;
#line 198 "slasq2.f"
	    s = z__[3] * (z__[2] / t);
#line 199 "slasq2.f"
	    if (s <= t) {
#line 200 "slasq2.f"
		s = z__[3] * (z__[2] / (t * (sqrt(s / t + 1.) + 1.)));
#line 201 "slasq2.f"
	    } else {
#line 202 "slasq2.f"
		s = z__[3] * (z__[2] / (t + sqrt(t) * sqrt(t + s)));
#line 203 "slasq2.f"
	    }
#line 204 "slasq2.f"
	    t = z__[1] + (s + z__[2]);
#line 205 "slasq2.f"
	    z__[3] *= z__[1] / t;
#line 206 "slasq2.f"
	    z__[1] = t;
#line 207 "slasq2.f"
	}
#line 208 "slasq2.f"
	z__[2] = z__[3];
#line 209 "slasq2.f"
	z__[6] = z__[2] + z__[1];
#line 210 "slasq2.f"
	return 0;
#line 211 "slasq2.f"
    }

/*     Check for negative data and compute sums of q's and e's. */

#line 215 "slasq2.f"
    z__[*n * 2] = 0.;
#line 216 "slasq2.f"
    emin = z__[2];
#line 217 "slasq2.f"
    qmax = 0.;
#line 218 "slasq2.f"
    zmax = 0.;
#line 219 "slasq2.f"
    d__ = 0.;
#line 220 "slasq2.f"
    e = 0.;

#line 222 "slasq2.f"
    i__1 = *n - 1 << 1;
#line 222 "slasq2.f"
    for (k = 1; k <= i__1; k += 2) {
#line 223 "slasq2.f"
	if (z__[k] < 0.) {
#line 224 "slasq2.f"
	    *info = -(k + 200);
#line 225 "slasq2.f"
	    xerbla_("SLASQ2", &c__2, (ftnlen)6);
#line 226 "slasq2.f"
	    return 0;
#line 227 "slasq2.f"
	} else if (z__[k + 1] < 0.) {
#line 228 "slasq2.f"
	    *info = -(k + 201);
#line 229 "slasq2.f"
	    xerbla_("SLASQ2", &c__2, (ftnlen)6);
#line 230 "slasq2.f"
	    return 0;
#line 231 "slasq2.f"
	}
#line 232 "slasq2.f"
	d__ += z__[k];
#line 233 "slasq2.f"
	e += z__[k + 1];
/* Computing MAX */
#line 234 "slasq2.f"
	d__1 = qmax, d__2 = z__[k];
#line 234 "slasq2.f"
	qmax = max(d__1,d__2);
/* Computing MIN */
#line 235 "slasq2.f"
	d__1 = emin, d__2 = z__[k + 1];
#line 235 "slasq2.f"
	emin = min(d__1,d__2);
/* Computing MAX */
#line 236 "slasq2.f"
	d__1 = max(qmax,zmax), d__2 = z__[k + 1];
#line 236 "slasq2.f"
	zmax = max(d__1,d__2);
#line 237 "slasq2.f"
/* L10: */
#line 237 "slasq2.f"
    }
#line 238 "slasq2.f"
    if (z__[(*n << 1) - 1] < 0.) {
#line 239 "slasq2.f"
	*info = -((*n << 1) + 199);
#line 240 "slasq2.f"
	xerbla_("SLASQ2", &c__2, (ftnlen)6);
#line 241 "slasq2.f"
	return 0;
#line 242 "slasq2.f"
    }
#line 243 "slasq2.f"
    d__ += z__[(*n << 1) - 1];
/* Computing MAX */
#line 244 "slasq2.f"
    d__1 = qmax, d__2 = z__[(*n << 1) - 1];
#line 244 "slasq2.f"
    qmax = max(d__1,d__2);
#line 245 "slasq2.f"
    zmax = max(qmax,zmax);

/*     Check for diagonality. */

#line 249 "slasq2.f"
    if (e == 0.) {
#line 250 "slasq2.f"
	i__1 = *n;
#line 250 "slasq2.f"
	for (k = 2; k <= i__1; ++k) {
#line 251 "slasq2.f"
	    z__[k] = z__[(k << 1) - 1];
#line 252 "slasq2.f"
/* L20: */
#line 252 "slasq2.f"
	}
#line 253 "slasq2.f"
	slasrt_("D", n, &z__[1], &iinfo, (ftnlen)1);
#line 254 "slasq2.f"
	z__[(*n << 1) - 1] = d__;
#line 255 "slasq2.f"
	return 0;
#line 256 "slasq2.f"
    }

#line 258 "slasq2.f"
    trace = d__ + e;

/*     Check for zero data. */

#line 262 "slasq2.f"
    if (trace == 0.) {
#line 263 "slasq2.f"
	z__[(*n << 1) - 1] = 0.;
#line 264 "slasq2.f"
	return 0;
#line 265 "slasq2.f"
    }

/*     Check whether the machine is IEEE conformable. */

/*     IEEE = ILAENV( 10, 'SLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 .AND. */
/*    $       ILAENV( 11, 'SLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 */

/*     [11/15/2008] The case IEEE=.TRUE. has a problem in single precision with */
/*     some the test matrices of type 16. The double precision code is fine. */

#line 275 "slasq2.f"
    ieee = FALSE_;

/*     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...). */

#line 279 "slasq2.f"
    for (k = *n << 1; k >= 2; k += -2) {
#line 280 "slasq2.f"
	z__[k * 2] = 0.;
#line 281 "slasq2.f"
	z__[(k << 1) - 1] = z__[k];
#line 282 "slasq2.f"
	z__[(k << 1) - 2] = 0.;
#line 283 "slasq2.f"
	z__[(k << 1) - 3] = z__[k - 1];
#line 284 "slasq2.f"
/* L30: */
#line 284 "slasq2.f"
    }

#line 286 "slasq2.f"
    i0 = 1;
#line 287 "slasq2.f"
    n0 = *n;

/*     Reverse the qd-array, if warranted. */

#line 291 "slasq2.f"
    if (z__[(i0 << 2) - 3] * 1.5 < z__[(n0 << 2) - 3]) {
#line 292 "slasq2.f"
	ipn4 = i0 + n0 << 2;
#line 293 "slasq2.f"
	i__1 = i0 + n0 - 1 << 1;
#line 293 "slasq2.f"
	for (i4 = i0 << 2; i4 <= i__1; i4 += 4) {
#line 294 "slasq2.f"
	    temp = z__[i4 - 3];
#line 295 "slasq2.f"
	    z__[i4 - 3] = z__[ipn4 - i4 - 3];
#line 296 "slasq2.f"
	    z__[ipn4 - i4 - 3] = temp;
#line 297 "slasq2.f"
	    temp = z__[i4 - 1];
#line 298 "slasq2.f"
	    z__[i4 - 1] = z__[ipn4 - i4 - 5];
#line 299 "slasq2.f"
	    z__[ipn4 - i4 - 5] = temp;
#line 300 "slasq2.f"
/* L40: */
#line 300 "slasq2.f"
	}
#line 301 "slasq2.f"
    }

/*     Initial split checking via dqd and Li's test. */

#line 305 "slasq2.f"
    pp = 0;

#line 307 "slasq2.f"
    for (k = 1; k <= 2; ++k) {

#line 309 "slasq2.f"
	d__ = z__[(n0 << 2) + pp - 3];
#line 310 "slasq2.f"
	i__1 = (i0 << 2) + pp;
#line 310 "slasq2.f"
	for (i4 = (n0 - 1 << 2) + pp; i4 >= i__1; i4 += -4) {
#line 311 "slasq2.f"
	    if (z__[i4 - 1] <= tol2 * d__) {
#line 312 "slasq2.f"
		z__[i4 - 1] = -0.;
#line 313 "slasq2.f"
		d__ = z__[i4 - 3];
#line 314 "slasq2.f"
	    } else {
#line 315 "slasq2.f"
		d__ = z__[i4 - 3] * (d__ / (d__ + z__[i4 - 1]));
#line 316 "slasq2.f"
	    }
#line 317 "slasq2.f"
/* L50: */
#line 317 "slasq2.f"
	}

/*        dqd maps Z to ZZ plus Li's test. */

#line 321 "slasq2.f"
	emin = z__[(i0 << 2) + pp + 1];
#line 322 "slasq2.f"
	d__ = z__[(i0 << 2) + pp - 3];
#line 323 "slasq2.f"
	i__1 = (n0 - 1 << 2) + pp;
#line 323 "slasq2.f"
	for (i4 = (i0 << 2) + pp; i4 <= i__1; i4 += 4) {
#line 324 "slasq2.f"
	    z__[i4 - (pp << 1) - 2] = d__ + z__[i4 - 1];
#line 325 "slasq2.f"
	    if (z__[i4 - 1] <= tol2 * d__) {
#line 326 "slasq2.f"
		z__[i4 - 1] = -0.;
#line 327 "slasq2.f"
		z__[i4 - (pp << 1) - 2] = d__;
#line 328 "slasq2.f"
		z__[i4 - (pp << 1)] = 0.;
#line 329 "slasq2.f"
		d__ = z__[i4 + 1];
#line 330 "slasq2.f"
	    } else if (safmin * z__[i4 + 1] < z__[i4 - (pp << 1) - 2] && 
		    safmin * z__[i4 - (pp << 1) - 2] < z__[i4 + 1]) {
#line 332 "slasq2.f"
		temp = z__[i4 + 1] / z__[i4 - (pp << 1) - 2];
#line 333 "slasq2.f"
		z__[i4 - (pp << 1)] = z__[i4 - 1] * temp;
#line 334 "slasq2.f"
		d__ *= temp;
#line 335 "slasq2.f"
	    } else {
#line 336 "slasq2.f"
		z__[i4 - (pp << 1)] = z__[i4 + 1] * (z__[i4 - 1] / z__[i4 - (
			pp << 1) - 2]);
#line 337 "slasq2.f"
		d__ = z__[i4 + 1] * (d__ / z__[i4 - (pp << 1) - 2]);
#line 338 "slasq2.f"
	    }
/* Computing MIN */
#line 339 "slasq2.f"
	    d__1 = emin, d__2 = z__[i4 - (pp << 1)];
#line 339 "slasq2.f"
	    emin = min(d__1,d__2);
#line 340 "slasq2.f"
/* L60: */
#line 340 "slasq2.f"
	}
#line 341 "slasq2.f"
	z__[(n0 << 2) - pp - 2] = d__;

/*        Now find qmax. */

#line 345 "slasq2.f"
	qmax = z__[(i0 << 2) - pp - 2];
#line 346 "slasq2.f"
	i__1 = (n0 << 2) - pp - 2;
#line 346 "slasq2.f"
	for (i4 = (i0 << 2) - pp + 2; i4 <= i__1; i4 += 4) {
/* Computing MAX */
#line 347 "slasq2.f"
	    d__1 = qmax, d__2 = z__[i4];
#line 347 "slasq2.f"
	    qmax = max(d__1,d__2);
#line 348 "slasq2.f"
/* L70: */
#line 348 "slasq2.f"
	}

/*        Prepare for the next iteration on K. */

#line 352 "slasq2.f"
	pp = 1 - pp;
#line 353 "slasq2.f"
/* L80: */
#line 353 "slasq2.f"
    }

/*     Initialise variables to pass to SLASQ3. */

#line 357 "slasq2.f"
    ttype = 0;
#line 358 "slasq2.f"
    dmin1 = 0.;
#line 359 "slasq2.f"
    dmin2 = 0.;
#line 360 "slasq2.f"
    dn = 0.;
#line 361 "slasq2.f"
    dn1 = 0.;
#line 362 "slasq2.f"
    dn2 = 0.;
#line 363 "slasq2.f"
    g = 0.;
#line 364 "slasq2.f"
    tau = 0.;

#line 366 "slasq2.f"
    iter = 2;
#line 367 "slasq2.f"
    nfail = 0;
#line 368 "slasq2.f"
    ndiv = n0 - i0 << 1;

#line 370 "slasq2.f"
    i__1 = *n + 1;
#line 370 "slasq2.f"
    for (iwhila = 1; iwhila <= i__1; ++iwhila) {
#line 371 "slasq2.f"
	if (n0 < 1) {
#line 371 "slasq2.f"
	    goto L170;
#line 371 "slasq2.f"
	}

/*        While array unfinished do */

/*        E(N0) holds the value of SIGMA when submatrix in I0:N0 */
/*        splits from the rest of the array, but is negated. */

#line 379 "slasq2.f"
	desig = 0.;
#line 380 "slasq2.f"
	if (n0 == *n) {
#line 381 "slasq2.f"
	    sigma = 0.;
#line 382 "slasq2.f"
	} else {
#line 383 "slasq2.f"
	    sigma = -z__[(n0 << 2) - 1];
#line 384 "slasq2.f"
	}
#line 385 "slasq2.f"
	if (sigma < 0.) {
#line 386 "slasq2.f"
	    *info = 1;
#line 387 "slasq2.f"
	    return 0;
#line 388 "slasq2.f"
	}

/*        Find last unreduced submatrix's top index I0, find QMAX and */
/*        EMIN. Find Gershgorin-type bound if Q's much greater than E's. */

#line 393 "slasq2.f"
	emax = 0.;
#line 394 "slasq2.f"
	if (n0 > i0) {
#line 395 "slasq2.f"
	    emin = (d__1 = z__[(n0 << 2) - 5], abs(d__1));
#line 396 "slasq2.f"
	} else {
#line 397 "slasq2.f"
	    emin = 0.;
#line 398 "slasq2.f"
	}
#line 399 "slasq2.f"
	qmin = z__[(n0 << 2) - 3];
#line 400 "slasq2.f"
	qmax = qmin;
#line 401 "slasq2.f"
	for (i4 = n0 << 2; i4 >= 8; i4 += -4) {
#line 402 "slasq2.f"
	    if (z__[i4 - 5] <= 0.) {
#line 402 "slasq2.f"
		goto L100;
#line 402 "slasq2.f"
	    }
#line 404 "slasq2.f"
	    if (qmin >= emax * 4.) {
/* Computing MIN */
#line 405 "slasq2.f"
		d__1 = qmin, d__2 = z__[i4 - 3];
#line 405 "slasq2.f"
		qmin = min(d__1,d__2);
/* Computing MAX */
#line 406 "slasq2.f"
		d__1 = emax, d__2 = z__[i4 - 5];
#line 406 "slasq2.f"
		emax = max(d__1,d__2);
#line 407 "slasq2.f"
	    }
/* Computing MAX */
#line 408 "slasq2.f"
	    d__1 = qmax, d__2 = z__[i4 - 7] + z__[i4 - 5];
#line 408 "slasq2.f"
	    qmax = max(d__1,d__2);
/* Computing MIN */
#line 409 "slasq2.f"
	    d__1 = emin, d__2 = z__[i4 - 5];
#line 409 "slasq2.f"
	    emin = min(d__1,d__2);
#line 410 "slasq2.f"
/* L90: */
#line 410 "slasq2.f"
	}
#line 411 "slasq2.f"
	i4 = 4;

#line 413 "slasq2.f"
L100:
#line 414 "slasq2.f"
	i0 = i4 / 4;
#line 415 "slasq2.f"
	pp = 0;

#line 417 "slasq2.f"
	if (n0 - i0 > 1) {
#line 418 "slasq2.f"
	    dee = z__[(i0 << 2) - 3];
#line 419 "slasq2.f"
	    deemin = dee;
#line 420 "slasq2.f"
	    kmin = i0;
#line 421 "slasq2.f"
	    i__2 = (n0 << 2) - 3;
#line 421 "slasq2.f"
	    for (i4 = (i0 << 2) + 1; i4 <= i__2; i4 += 4) {
#line 422 "slasq2.f"
		dee = z__[i4] * (dee / (dee + z__[i4 - 2]));
#line 423 "slasq2.f"
		if (dee <= deemin) {
#line 424 "slasq2.f"
		    deemin = dee;
#line 425 "slasq2.f"
		    kmin = (i4 + 3) / 4;
#line 426 "slasq2.f"
		}
#line 427 "slasq2.f"
/* L110: */
#line 427 "slasq2.f"
	    }
#line 428 "slasq2.f"
	    if (kmin - i0 << 1 < n0 - kmin && deemin <= z__[(n0 << 2) - 3] * 
		    .5) {
#line 430 "slasq2.f"
		ipn4 = i0 + n0 << 2;
#line 431 "slasq2.f"
		pp = 2;
#line 432 "slasq2.f"
		i__2 = i0 + n0 - 1 << 1;
#line 432 "slasq2.f"
		for (i4 = i0 << 2; i4 <= i__2; i4 += 4) {
#line 433 "slasq2.f"
		    temp = z__[i4 - 3];
#line 434 "slasq2.f"
		    z__[i4 - 3] = z__[ipn4 - i4 - 3];
#line 435 "slasq2.f"
		    z__[ipn4 - i4 - 3] = temp;
#line 436 "slasq2.f"
		    temp = z__[i4 - 2];
#line 437 "slasq2.f"
		    z__[i4 - 2] = z__[ipn4 - i4 - 2];
#line 438 "slasq2.f"
		    z__[ipn4 - i4 - 2] = temp;
#line 439 "slasq2.f"
		    temp = z__[i4 - 1];
#line 440 "slasq2.f"
		    z__[i4 - 1] = z__[ipn4 - i4 - 5];
#line 441 "slasq2.f"
		    z__[ipn4 - i4 - 5] = temp;
#line 442 "slasq2.f"
		    temp = z__[i4];
#line 443 "slasq2.f"
		    z__[i4] = z__[ipn4 - i4 - 4];
#line 444 "slasq2.f"
		    z__[ipn4 - i4 - 4] = temp;
#line 445 "slasq2.f"
/* L120: */
#line 445 "slasq2.f"
		}
#line 446 "slasq2.f"
	    }
#line 447 "slasq2.f"
	}

/*        Put -(initial shift) into DMIN. */

/* Computing MAX */
#line 451 "slasq2.f"
	d__1 = 0., d__2 = qmin - sqrt(qmin) * 2. * sqrt(emax);
#line 451 "slasq2.f"
	dmin__ = -max(d__1,d__2);

/*        Now I0:N0 is unreduced. */
/*        PP = 0 for ping, PP = 1 for pong. */
/*        PP = 2 indicates that flipping was applied to the Z array and */
/*               and that the tests for deflation upon entry in SLASQ3 */
/*               should not be performed. */

#line 459 "slasq2.f"
	nbig = (n0 - i0 + 1) * 100;
#line 460 "slasq2.f"
	i__2 = nbig;
#line 460 "slasq2.f"
	for (iwhilb = 1; iwhilb <= i__2; ++iwhilb) {
#line 461 "slasq2.f"
	    if (i0 > n0) {
#line 461 "slasq2.f"
		goto L150;
#line 461 "slasq2.f"
	    }

/*           While submatrix unfinished take a good dqds step. */

#line 466 "slasq2.f"
	    slasq3_(&i0, &n0, &z__[1], &pp, &dmin__, &sigma, &desig, &qmax, &
		    nfail, &iter, &ndiv, &ieee, &ttype, &dmin1, &dmin2, &dn, &
		    dn1, &dn2, &g, &tau);

#line 470 "slasq2.f"
	    pp = 1 - pp;

/*           When EMIN is very small check for splits. */

#line 474 "slasq2.f"
	    if (pp == 0 && n0 - i0 >= 3) {
#line 475 "slasq2.f"
		if (z__[n0 * 4] <= tol2 * qmax || z__[(n0 << 2) - 1] <= tol2 *
			 sigma) {
#line 477 "slasq2.f"
		    splt = i0 - 1;
#line 478 "slasq2.f"
		    qmax = z__[(i0 << 2) - 3];
#line 479 "slasq2.f"
		    emin = z__[(i0 << 2) - 1];
#line 480 "slasq2.f"
		    oldemn = z__[i0 * 4];
#line 481 "slasq2.f"
		    i__3 = n0 - 3 << 2;
#line 481 "slasq2.f"
		    for (i4 = i0 << 2; i4 <= i__3; i4 += 4) {
#line 482 "slasq2.f"
			if (z__[i4] <= tol2 * z__[i4 - 3] || z__[i4 - 1] <= 
				tol2 * sigma) {
#line 484 "slasq2.f"
			    z__[i4 - 1] = -sigma;
#line 485 "slasq2.f"
			    splt = i4 / 4;
#line 486 "slasq2.f"
			    qmax = 0.;
#line 487 "slasq2.f"
			    emin = z__[i4 + 3];
#line 488 "slasq2.f"
			    oldemn = z__[i4 + 4];
#line 489 "slasq2.f"
			} else {
/* Computing MAX */
#line 490 "slasq2.f"
			    d__1 = qmax, d__2 = z__[i4 + 1];
#line 490 "slasq2.f"
			    qmax = max(d__1,d__2);
/* Computing MIN */
#line 491 "slasq2.f"
			    d__1 = emin, d__2 = z__[i4 - 1];
#line 491 "slasq2.f"
			    emin = min(d__1,d__2);
/* Computing MIN */
#line 492 "slasq2.f"
			    d__1 = oldemn, d__2 = z__[i4];
#line 492 "slasq2.f"
			    oldemn = min(d__1,d__2);
#line 493 "slasq2.f"
			}
#line 494 "slasq2.f"
/* L130: */
#line 494 "slasq2.f"
		    }
#line 495 "slasq2.f"
		    z__[(n0 << 2) - 1] = emin;
#line 496 "slasq2.f"
		    z__[n0 * 4] = oldemn;
#line 497 "slasq2.f"
		    i0 = splt + 1;
#line 498 "slasq2.f"
		}
#line 499 "slasq2.f"
	    }

#line 501 "slasq2.f"
/* L140: */
#line 501 "slasq2.f"
	}

#line 503 "slasq2.f"
	*info = 2;

/*        Maximum number of iterations exceeded, restore the shift */
/*        SIGMA and place the new d's and e's in a qd array. */
/*        This might need to be done for several blocks */

#line 509 "slasq2.f"
	i1 = i0;
#line 510 "slasq2.f"
	n1 = n0;
#line 511 "slasq2.f"
L145:
#line 512 "slasq2.f"
	tempq = z__[(i0 << 2) - 3];
#line 513 "slasq2.f"
	z__[(i0 << 2) - 3] += sigma;
#line 514 "slasq2.f"
	i__2 = n0;
#line 514 "slasq2.f"
	for (k = i0 + 1; k <= i__2; ++k) {
#line 515 "slasq2.f"
	    tempe = z__[(k << 2) - 5];
#line 516 "slasq2.f"
	    z__[(k << 2) - 5] *= tempq / z__[(k << 2) - 7];
#line 517 "slasq2.f"
	    tempq = z__[(k << 2) - 3];
#line 518 "slasq2.f"
	    z__[(k << 2) - 3] = z__[(k << 2) - 3] + sigma + tempe - z__[(k << 
		    2) - 5];
#line 519 "slasq2.f"
	}

/*        Prepare to do this on the previous block if there is one */

#line 523 "slasq2.f"
	if (i1 > 1) {
#line 524 "slasq2.f"
	    n1 = i1 - 1;
#line 525 "slasq2.f"
	    while(i1 >= 2 && z__[(i1 << 2) - 5] >= 0.) {
#line 526 "slasq2.f"
		--i1;
#line 527 "slasq2.f"
	    }
#line 528 "slasq2.f"
	    if (i1 >= 1) {
#line 529 "slasq2.f"
		sigma = -z__[(n1 << 2) - 1];
#line 530 "slasq2.f"
		goto L145;
#line 531 "slasq2.f"
	    }
#line 532 "slasq2.f"
	}
#line 534 "slasq2.f"
	i__2 = *n;
#line 534 "slasq2.f"
	for (k = 1; k <= i__2; ++k) {
#line 535 "slasq2.f"
	    z__[(k << 1) - 1] = z__[(k << 2) - 3];

/*        Only the block 1..N0 is unfinished.  The rest of the e's */
/*        must be essentially zero, although sometimes other data */
/*        has been stored in them. */

#line 541 "slasq2.f"
	    if (k < n0) {
#line 542 "slasq2.f"
		z__[k * 2] = z__[(k << 2) - 1];
#line 543 "slasq2.f"
	    } else {
#line 544 "slasq2.f"
		z__[k * 2] = 0.;
#line 545 "slasq2.f"
	    }
#line 546 "slasq2.f"
	}
#line 547 "slasq2.f"
	return 0;

/*        end IWHILB */

#line 551 "slasq2.f"
L150:

#line 553 "slasq2.f"
/* L160: */
#line 553 "slasq2.f"
	;
#line 553 "slasq2.f"
    }

#line 555 "slasq2.f"
    *info = 3;
#line 556 "slasq2.f"
    return 0;

/*     end IWHILA */

#line 560 "slasq2.f"
L170:

/*     Move q's to the front. */

#line 564 "slasq2.f"
    i__1 = *n;
#line 564 "slasq2.f"
    for (k = 2; k <= i__1; ++k) {
#line 565 "slasq2.f"
	z__[k] = z__[(k << 2) - 3];
#line 566 "slasq2.f"
/* L180: */
#line 566 "slasq2.f"
    }

/*     Sort and compute sum of eigenvalues. */

#line 570 "slasq2.f"
    slasrt_("D", n, &z__[1], &iinfo, (ftnlen)1);

#line 572 "slasq2.f"
    e = 0.;
#line 573 "slasq2.f"
    for (k = *n; k >= 1; --k) {
#line 574 "slasq2.f"
	e += z__[k];
#line 575 "slasq2.f"
/* L190: */
#line 575 "slasq2.f"
    }

/*     Store trace, sum(eigenvalues) and information on performance. */

#line 579 "slasq2.f"
    z__[(*n << 1) + 1] = trace;
#line 580 "slasq2.f"
    z__[(*n << 1) + 2] = e;
#line 581 "slasq2.f"
    z__[(*n << 1) + 3] = (doublereal) iter;
/* Computing 2nd power */
#line 582 "slasq2.f"
    i__1 = *n;
#line 582 "slasq2.f"
    z__[(*n << 1) + 4] = (doublereal) ndiv / (doublereal) (i__1 * i__1);
#line 583 "slasq2.f"
    z__[(*n << 1) + 5] = nfail * 100. / (doublereal) iter;
#line 584 "slasq2.f"
    return 0;

/*     End of SLASQ2 */

} /* slasq2_ */

