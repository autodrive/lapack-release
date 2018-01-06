#line 1 "slasy2.f"
/* slasy2.f -- translated by f2c (version 20100827).
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

#line 1 "slasy2.f"
/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__16 = 16;
static integer c__0 = 0;

/* > \brief \b SLASY2 solves the Sylvester matrix equation where the matrices are of order 1 or 2. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASY2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasy2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasy2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasy2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, */
/*                          LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            LTRANL, LTRANR */
/*       INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2 */
/*       REAL               SCALE, XNORM */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in */
/* > */
/* >        op(TL)*X + ISGN*X*op(TR) = SCALE*B, */
/* > */
/* > where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or */
/* > -1.  op(T) = T or T**T, where T**T denotes the transpose of T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] LTRANL */
/* > \verbatim */
/* >          LTRANL is LOGICAL */
/* >          On entry, LTRANL specifies the op(TL): */
/* >             = .FALSE., op(TL) = TL, */
/* >             = .TRUE., op(TL) = TL**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LTRANR */
/* > \verbatim */
/* >          LTRANR is LOGICAL */
/* >          On entry, LTRANR specifies the op(TR): */
/* >            = .FALSE., op(TR) = TR, */
/* >            = .TRUE., op(TR) = TR**T. */
/* > \endverbatim */
/* > */
/* > \param[in] ISGN */
/* > \verbatim */
/* >          ISGN is INTEGER */
/* >          On entry, ISGN specifies the sign of the equation */
/* >          as described before. ISGN may only be 1 or -1. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* >          On entry, N1 specifies the order of matrix TL. */
/* >          N1 may only be 0, 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* >          N2 is INTEGER */
/* >          On entry, N2 specifies the order of matrix TR. */
/* >          N2 may only be 0, 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] TL */
/* > \verbatim */
/* >          TL is REAL array, dimension (LDTL,2) */
/* >          On entry, TL contains an N1 by N1 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDTL */
/* > \verbatim */
/* >          LDTL is INTEGER */
/* >          The leading dimension of the matrix TL. LDTL >= max(1,N1). */
/* > \endverbatim */
/* > */
/* > \param[in] TR */
/* > \verbatim */
/* >          TR is REAL array, dimension (LDTR,2) */
/* >          On entry, TR contains an N2 by N2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDTR */
/* > \verbatim */
/* >          LDTR is INTEGER */
/* >          The leading dimension of the matrix TR. LDTR >= max(1,N2). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,2) */
/* >          On entry, the N1 by N2 matrix B contains the right-hand */
/* >          side of the equation. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the matrix B. LDB >= max(1,N1). */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >          On exit, SCALE contains the scale factor. SCALE is chosen */
/* >          less than or equal to 1 to prevent the solution overflowing. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (LDX,2) */
/* >          On exit, X contains the N1 by N2 solution. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the matrix X. LDX >= max(1,N1). */
/* > \endverbatim */
/* > */
/* > \param[out] XNORM */
/* > \verbatim */
/* >          XNORM is REAL */
/* >          On exit, XNORM is the infinity-norm of the solution. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          On exit, INFO is set to */
/* >             0: successful exit. */
/* >             1: TL and TR have too close eigenvalues, so TL or */
/* >                TR is perturbed to get a nonsingular equation. */
/* >          NOTE: In the interests of speed, this routine does not */
/* >                check the inputs for errors. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realSYauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slasy2_(logical *ltranl, logical *ltranr, integer *isgn, 
	integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
	tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale, 
	doublereal *x, integer *ldx, doublereal *xnorm, integer *info)
{
    /* Initialized data */

    static integer locu12[4] = { 3,4,1,2 };
    static integer locl21[4] = { 2,1,4,3 };
    static integer locu22[4] = { 4,3,2,1 };
    static logical xswpiv[4] = { FALSE_,FALSE_,TRUE_,TRUE_ };
    static logical bswpiv[4] = { FALSE_,TRUE_,FALSE_,TRUE_ };

    /* System generated locals */
    integer b_dim1, b_offset, tl_dim1, tl_offset, tr_dim1, tr_offset, x_dim1, 
	    x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Local variables */
    static integer i__, j, k;
    static doublereal x2[2], l21, u11, u12;
    static integer ip, jp;
    static doublereal u22, t16[16]	/* was [4][4] */, gam, bet, eps, sgn, 
	    tmp[4], tau1, btmp[4], smin;
    static integer ipiv;
    static doublereal temp;
    static integer jpiv[4];
    static doublereal xmax;
    static integer ipsv, jpsv;
    static logical bswap;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical xswap;
    extern doublereal slamch_(char *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
#line 224 "slasy2.f"
    /* Parameter adjustments */
#line 224 "slasy2.f"
    tl_dim1 = *ldtl;
#line 224 "slasy2.f"
    tl_offset = 1 + tl_dim1;
#line 224 "slasy2.f"
    tl -= tl_offset;
#line 224 "slasy2.f"
    tr_dim1 = *ldtr;
#line 224 "slasy2.f"
    tr_offset = 1 + tr_dim1;
#line 224 "slasy2.f"
    tr -= tr_offset;
#line 224 "slasy2.f"
    b_dim1 = *ldb;
#line 224 "slasy2.f"
    b_offset = 1 + b_dim1;
#line 224 "slasy2.f"
    b -= b_offset;
#line 224 "slasy2.f"
    x_dim1 = *ldx;
#line 224 "slasy2.f"
    x_offset = 1 + x_dim1;
#line 224 "slasy2.f"
    x -= x_offset;
#line 224 "slasy2.f"

#line 224 "slasy2.f"
    /* Function Body */
/*     .. */
/*     .. Executable Statements .. */

/*     Do not check the input parameters for errors */

#line 233 "slasy2.f"
    *info = 0;

/*     Quick return if possible */

#line 237 "slasy2.f"
    if (*n1 == 0 || *n2 == 0) {
#line 237 "slasy2.f"
	return 0;
#line 237 "slasy2.f"
    }

/*     Set constants to control overflow */

#line 242 "slasy2.f"
    eps = slamch_("P", (ftnlen)1);
#line 243 "slasy2.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 244 "slasy2.f"
    sgn = (doublereal) (*isgn);

#line 246 "slasy2.f"
    k = *n1 + *n1 + *n2 - 2;
#line 247 "slasy2.f"
    switch (k) {
#line 247 "slasy2.f"
	case 1:  goto L10;
#line 247 "slasy2.f"
	case 2:  goto L20;
#line 247 "slasy2.f"
	case 3:  goto L30;
#line 247 "slasy2.f"
	case 4:  goto L50;
#line 247 "slasy2.f"
    }

/*     1 by 1: TL11*X + SGN*X*TR11 = B11 */

#line 251 "slasy2.f"
L10:
#line 252 "slasy2.f"
    tau1 = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
#line 253 "slasy2.f"
    bet = abs(tau1);
#line 254 "slasy2.f"
    if (bet <= smlnum) {
#line 255 "slasy2.f"
	tau1 = smlnum;
#line 256 "slasy2.f"
	bet = smlnum;
#line 257 "slasy2.f"
	*info = 1;
#line 258 "slasy2.f"
    }

#line 260 "slasy2.f"
    *scale = 1.;
#line 261 "slasy2.f"
    gam = (d__1 = b[b_dim1 + 1], abs(d__1));
#line 262 "slasy2.f"
    if (smlnum * gam > bet) {
#line 262 "slasy2.f"
	*scale = 1. / gam;
#line 262 "slasy2.f"
    }

#line 265 "slasy2.f"
    x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / tau1;
#line 266 "slasy2.f"
    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1));
#line 267 "slasy2.f"
    return 0;

/*     1 by 2: */
/*     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12] */
/*                                       [TR21 TR22] */

#line 273 "slasy2.f"
L20:

/* Computing MAX */
/* Computing MAX */
#line 275 "slasy2.f"
    d__7 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__8 = (d__2 = tr[tr_dim1 + 1]
	    , abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = tr[(tr_dim1 <<
	     1) + 1], abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = tr[
	    tr_dim1 + 2], abs(d__4)), d__7 = max(d__7,d__8), d__8 = (d__5 = 
	    tr[(tr_dim1 << 1) + 2], abs(d__5));
#line 275 "slasy2.f"
    d__6 = eps * max(d__7,d__8);
#line 275 "slasy2.f"
    smin = max(d__6,smlnum);
#line 278 "slasy2.f"
    tmp[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
#line 279 "slasy2.f"
    tmp[3] = tl[tl_dim1 + 1] + sgn * tr[(tr_dim1 << 1) + 2];
#line 280 "slasy2.f"
    if (*ltranr) {
#line 281 "slasy2.f"
	tmp[1] = sgn * tr[tr_dim1 + 2];
#line 282 "slasy2.f"
	tmp[2] = sgn * tr[(tr_dim1 << 1) + 1];
#line 283 "slasy2.f"
    } else {
#line 284 "slasy2.f"
	tmp[1] = sgn * tr[(tr_dim1 << 1) + 1];
#line 285 "slasy2.f"
	tmp[2] = sgn * tr[tr_dim1 + 2];
#line 286 "slasy2.f"
    }
#line 287 "slasy2.f"
    btmp[0] = b[b_dim1 + 1];
#line 288 "slasy2.f"
    btmp[1] = b[(b_dim1 << 1) + 1];
#line 289 "slasy2.f"
    goto L40;

/*     2 by 1: */
/*          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11] */
/*            [TL21 TL22] [X21]         [X21]         [B21] */

#line 295 "slasy2.f"
L30:
/* Computing MAX */
/* Computing MAX */
#line 296 "slasy2.f"
    d__7 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__8 = (d__2 = tl[tl_dim1 + 1]
	    , abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = tl[(tl_dim1 <<
	     1) + 1], abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = tl[
	    tl_dim1 + 2], abs(d__4)), d__7 = max(d__7,d__8), d__8 = (d__5 = 
	    tl[(tl_dim1 << 1) + 2], abs(d__5));
#line 296 "slasy2.f"
    d__6 = eps * max(d__7,d__8);
#line 296 "slasy2.f"
    smin = max(d__6,smlnum);
#line 299 "slasy2.f"
    tmp[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
#line 300 "slasy2.f"
    tmp[3] = tl[(tl_dim1 << 1) + 2] + sgn * tr[tr_dim1 + 1];
#line 301 "slasy2.f"
    if (*ltranl) {
#line 302 "slasy2.f"
	tmp[1] = tl[(tl_dim1 << 1) + 1];
#line 303 "slasy2.f"
	tmp[2] = tl[tl_dim1 + 2];
#line 304 "slasy2.f"
    } else {
#line 305 "slasy2.f"
	tmp[1] = tl[tl_dim1 + 2];
#line 306 "slasy2.f"
	tmp[2] = tl[(tl_dim1 << 1) + 1];
#line 307 "slasy2.f"
    }
#line 308 "slasy2.f"
    btmp[0] = b[b_dim1 + 1];
#line 309 "slasy2.f"
    btmp[1] = b[b_dim1 + 2];
#line 310 "slasy2.f"
L40:

/*     Solve 2 by 2 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 315 "slasy2.f"
    ipiv = isamax_(&c__4, tmp, &c__1);
#line 316 "slasy2.f"
    u11 = tmp[ipiv - 1];
#line 317 "slasy2.f"
    if (abs(u11) <= smin) {
#line 318 "slasy2.f"
	*info = 1;
#line 319 "slasy2.f"
	u11 = smin;
#line 320 "slasy2.f"
    }
#line 321 "slasy2.f"
    u12 = tmp[locu12[ipiv - 1] - 1];
#line 322 "slasy2.f"
    l21 = tmp[locl21[ipiv - 1] - 1] / u11;
#line 323 "slasy2.f"
    u22 = tmp[locu22[ipiv - 1] - 1] - u12 * l21;
#line 324 "slasy2.f"
    xswap = xswpiv[ipiv - 1];
#line 325 "slasy2.f"
    bswap = bswpiv[ipiv - 1];
#line 326 "slasy2.f"
    if (abs(u22) <= smin) {
#line 327 "slasy2.f"
	*info = 1;
#line 328 "slasy2.f"
	u22 = smin;
#line 329 "slasy2.f"
    }
#line 330 "slasy2.f"
    if (bswap) {
#line 331 "slasy2.f"
	temp = btmp[1];
#line 332 "slasy2.f"
	btmp[1] = btmp[0] - l21 * temp;
#line 333 "slasy2.f"
	btmp[0] = temp;
#line 334 "slasy2.f"
    } else {
#line 335 "slasy2.f"
	btmp[1] -= l21 * btmp[0];
#line 336 "slasy2.f"
    }
#line 337 "slasy2.f"
    *scale = 1.;
#line 338 "slasy2.f"
    if (smlnum * 2. * abs(btmp[1]) > abs(u22) || smlnum * 2. * abs(btmp[0]) > 
	    abs(u11)) {
/* Computing MAX */
#line 340 "slasy2.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]);
#line 340 "slasy2.f"
	*scale = .5 / max(d__1,d__2);
#line 341 "slasy2.f"
	btmp[0] *= *scale;
#line 342 "slasy2.f"
	btmp[1] *= *scale;
#line 343 "slasy2.f"
    }
#line 344 "slasy2.f"
    x2[1] = btmp[1] / u22;
#line 345 "slasy2.f"
    x2[0] = btmp[0] / u11 - u12 / u11 * x2[1];
#line 346 "slasy2.f"
    if (xswap) {
#line 347 "slasy2.f"
	temp = x2[1];
#line 348 "slasy2.f"
	x2[1] = x2[0];
#line 349 "slasy2.f"
	x2[0] = temp;
#line 350 "slasy2.f"
    }
#line 351 "slasy2.f"
    x[x_dim1 + 1] = x2[0];
#line 352 "slasy2.f"
    if (*n1 == 1) {
#line 353 "slasy2.f"
	x[(x_dim1 << 1) + 1] = x2[1];
#line 354 "slasy2.f"
	*xnorm = (d__1 = x[x_dim1 + 1], abs(d__1)) + (d__2 = x[(x_dim1 << 1) 
		+ 1], abs(d__2));
#line 355 "slasy2.f"
    } else {
#line 356 "slasy2.f"
	x[x_dim1 + 2] = x2[1];
/* Computing MAX */
#line 357 "slasy2.f"
	d__3 = (d__1 = x[x_dim1 + 1], abs(d__1)), d__4 = (d__2 = x[x_dim1 + 2]
		, abs(d__2));
#line 357 "slasy2.f"
	*xnorm = max(d__3,d__4);
#line 358 "slasy2.f"
    }
#line 359 "slasy2.f"
    return 0;

/*     2 by 2: */
/*     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12] */
/*       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22] */

/*     Solve equivalent 4 by 4 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 368 "slasy2.f"
L50:
/* Computing MAX */
#line 369 "slasy2.f"
    d__5 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__6 = (d__2 = tr[(tr_dim1 << 
	    1) + 1], abs(d__2)), d__5 = max(d__5,d__6), d__6 = (d__3 = tr[
	    tr_dim1 + 2], abs(d__3)), d__5 = max(d__5,d__6), d__6 = (d__4 = 
	    tr[(tr_dim1 << 1) + 2], abs(d__4));
#line 369 "slasy2.f"
    smin = max(d__5,d__6);
/* Computing MAX */
#line 371 "slasy2.f"
    d__5 = smin, d__6 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__5 = max(d__5,
	    d__6), d__6 = (d__2 = tl[(tl_dim1 << 1) + 1], abs(d__2)), d__5 = 
	    max(d__5,d__6), d__6 = (d__3 = tl[tl_dim1 + 2], abs(d__3)), d__5 =
	     max(d__5,d__6), d__6 = (d__4 = tl[(tl_dim1 << 1) + 2], abs(d__4))
	    ;
#line 371 "slasy2.f"
    smin = max(d__5,d__6);
/* Computing MAX */
#line 373 "slasy2.f"
    d__1 = eps * smin;
#line 373 "slasy2.f"
    smin = max(d__1,smlnum);
#line 374 "slasy2.f"
    btmp[0] = 0.;
#line 375 "slasy2.f"
    scopy_(&c__16, btmp, &c__0, t16, &c__1);
#line 376 "slasy2.f"
    t16[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
#line 377 "slasy2.f"
    t16[5] = tl[(tl_dim1 << 1) + 2] + sgn * tr[tr_dim1 + 1];
#line 378 "slasy2.f"
    t16[10] = tl[tl_dim1 + 1] + sgn * tr[(tr_dim1 << 1) + 2];
#line 379 "slasy2.f"
    t16[15] = tl[(tl_dim1 << 1) + 2] + sgn * tr[(tr_dim1 << 1) + 2];
#line 380 "slasy2.f"
    if (*ltranl) {
#line 381 "slasy2.f"
	t16[4] = tl[tl_dim1 + 2];
#line 382 "slasy2.f"
	t16[1] = tl[(tl_dim1 << 1) + 1];
#line 383 "slasy2.f"
	t16[14] = tl[tl_dim1 + 2];
#line 384 "slasy2.f"
	t16[11] = tl[(tl_dim1 << 1) + 1];
#line 385 "slasy2.f"
    } else {
#line 386 "slasy2.f"
	t16[4] = tl[(tl_dim1 << 1) + 1];
#line 387 "slasy2.f"
	t16[1] = tl[tl_dim1 + 2];
#line 388 "slasy2.f"
	t16[14] = tl[(tl_dim1 << 1) + 1];
#line 389 "slasy2.f"
	t16[11] = tl[tl_dim1 + 2];
#line 390 "slasy2.f"
    }
#line 391 "slasy2.f"
    if (*ltranr) {
#line 392 "slasy2.f"
	t16[8] = sgn * tr[(tr_dim1 << 1) + 1];
#line 393 "slasy2.f"
	t16[13] = sgn * tr[(tr_dim1 << 1) + 1];
#line 394 "slasy2.f"
	t16[2] = sgn * tr[tr_dim1 + 2];
#line 395 "slasy2.f"
	t16[7] = sgn * tr[tr_dim1 + 2];
#line 396 "slasy2.f"
    } else {
#line 397 "slasy2.f"
	t16[8] = sgn * tr[tr_dim1 + 2];
#line 398 "slasy2.f"
	t16[13] = sgn * tr[tr_dim1 + 2];
#line 399 "slasy2.f"
	t16[2] = sgn * tr[(tr_dim1 << 1) + 1];
#line 400 "slasy2.f"
	t16[7] = sgn * tr[(tr_dim1 << 1) + 1];
#line 401 "slasy2.f"
    }
#line 402 "slasy2.f"
    btmp[0] = b[b_dim1 + 1];
#line 403 "slasy2.f"
    btmp[1] = b[b_dim1 + 2];
#line 404 "slasy2.f"
    btmp[2] = b[(b_dim1 << 1) + 1];
#line 405 "slasy2.f"
    btmp[3] = b[(b_dim1 << 1) + 2];

/*     Perform elimination */

#line 409 "slasy2.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 410 "slasy2.f"
	xmax = 0.;
#line 411 "slasy2.f"
	for (ip = i__; ip <= 4; ++ip) {
#line 412 "slasy2.f"
	    for (jp = i__; jp <= 4; ++jp) {
#line 413 "slasy2.f"
		if ((d__1 = t16[ip + (jp << 2) - 5], abs(d__1)) >= xmax) {
#line 414 "slasy2.f"
		    xmax = (d__1 = t16[ip + (jp << 2) - 5], abs(d__1));
#line 415 "slasy2.f"
		    ipsv = ip;
#line 416 "slasy2.f"
		    jpsv = jp;
#line 417 "slasy2.f"
		}
#line 418 "slasy2.f"
/* L60: */
#line 418 "slasy2.f"
	    }
#line 419 "slasy2.f"
/* L70: */
#line 419 "slasy2.f"
	}
#line 420 "slasy2.f"
	if (ipsv != i__) {
#line 421 "slasy2.f"
	    sswap_(&c__4, &t16[ipsv - 1], &c__4, &t16[i__ - 1], &c__4);
#line 422 "slasy2.f"
	    temp = btmp[i__ - 1];
#line 423 "slasy2.f"
	    btmp[i__ - 1] = btmp[ipsv - 1];
#line 424 "slasy2.f"
	    btmp[ipsv - 1] = temp;
#line 425 "slasy2.f"
	}
#line 426 "slasy2.f"
	if (jpsv != i__) {
#line 426 "slasy2.f"
	    sswap_(&c__4, &t16[(jpsv << 2) - 4], &c__1, &t16[(i__ << 2) - 4], 
		    &c__1);
#line 426 "slasy2.f"
	}
#line 428 "slasy2.f"
	jpiv[i__ - 1] = jpsv;
#line 429 "slasy2.f"
	if ((d__1 = t16[i__ + (i__ << 2) - 5], abs(d__1)) < smin) {
#line 430 "slasy2.f"
	    *info = 1;
#line 431 "slasy2.f"
	    t16[i__ + (i__ << 2) - 5] = smin;
#line 432 "slasy2.f"
	}
#line 433 "slasy2.f"
	for (j = i__ + 1; j <= 4; ++j) {
#line 434 "slasy2.f"
	    t16[j + (i__ << 2) - 5] /= t16[i__ + (i__ << 2) - 5];
#line 435 "slasy2.f"
	    btmp[j - 1] -= t16[j + (i__ << 2) - 5] * btmp[i__ - 1];
#line 436 "slasy2.f"
	    for (k = i__ + 1; k <= 4; ++k) {
#line 437 "slasy2.f"
		t16[j + (k << 2) - 5] -= t16[j + (i__ << 2) - 5] * t16[i__ + (
			k << 2) - 5];
#line 438 "slasy2.f"
/* L80: */
#line 438 "slasy2.f"
	    }
#line 439 "slasy2.f"
/* L90: */
#line 439 "slasy2.f"
	}
#line 440 "slasy2.f"
/* L100: */
#line 440 "slasy2.f"
    }
#line 441 "slasy2.f"
    if (abs(t16[15]) < smin) {
#line 441 "slasy2.f"
	t16[15] = smin;
#line 441 "slasy2.f"
    }
#line 443 "slasy2.f"
    *scale = 1.;
#line 444 "slasy2.f"
    if (smlnum * 8. * abs(btmp[0]) > abs(t16[0]) || smlnum * 8. * abs(btmp[1])
	     > abs(t16[5]) || smlnum * 8. * abs(btmp[2]) > abs(t16[10]) || 
	    smlnum * 8. * abs(btmp[3]) > abs(t16[15])) {
/* Computing MAX */
#line 448 "slasy2.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]), d__1 = max(d__1,d__2), d__2 
		= abs(btmp[2]), d__1 = max(d__1,d__2), d__2 = abs(btmp[3]);
#line 448 "slasy2.f"
	*scale = .125 / max(d__1,d__2);
#line 450 "slasy2.f"
	btmp[0] *= *scale;
#line 451 "slasy2.f"
	btmp[1] *= *scale;
#line 452 "slasy2.f"
	btmp[2] *= *scale;
#line 453 "slasy2.f"
	btmp[3] *= *scale;
#line 454 "slasy2.f"
    }
#line 455 "slasy2.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 456 "slasy2.f"
	k = 5 - i__;
#line 457 "slasy2.f"
	temp = 1. / t16[k + (k << 2) - 5];
#line 458 "slasy2.f"
	tmp[k - 1] = btmp[k - 1] * temp;
#line 459 "slasy2.f"
	for (j = k + 1; j <= 4; ++j) {
#line 460 "slasy2.f"
	    tmp[k - 1] -= temp * t16[k + (j << 2) - 5] * tmp[j - 1];
#line 461 "slasy2.f"
/* L110: */
#line 461 "slasy2.f"
	}
#line 462 "slasy2.f"
/* L120: */
#line 462 "slasy2.f"
    }
#line 463 "slasy2.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 464 "slasy2.f"
	if (jpiv[4 - i__ - 1] != 4 - i__) {
#line 465 "slasy2.f"
	    temp = tmp[4 - i__ - 1];
#line 466 "slasy2.f"
	    tmp[4 - i__ - 1] = tmp[jpiv[4 - i__ - 1] - 1];
#line 467 "slasy2.f"
	    tmp[jpiv[4 - i__ - 1] - 1] = temp;
#line 468 "slasy2.f"
	}
#line 469 "slasy2.f"
/* L130: */
#line 469 "slasy2.f"
    }
#line 470 "slasy2.f"
    x[x_dim1 + 1] = tmp[0];
#line 471 "slasy2.f"
    x[x_dim1 + 2] = tmp[1];
#line 472 "slasy2.f"
    x[(x_dim1 << 1) + 1] = tmp[2];
#line 473 "slasy2.f"
    x[(x_dim1 << 1) + 2] = tmp[3];
/* Computing MAX */
#line 474 "slasy2.f"
    d__1 = abs(tmp[0]) + abs(tmp[2]), d__2 = abs(tmp[1]) + abs(tmp[3]);
#line 474 "slasy2.f"
    *xnorm = max(d__1,d__2);
#line 476 "slasy2.f"
    return 0;

/*     End of SLASY2 */

} /* slasy2_ */

