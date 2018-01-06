#line 1 "dlasy2.f"
/* dlasy2.f -- translated by f2c (version 20100827).
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

#line 1 "dlasy2.f"
/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__16 = 16;
static integer c__0 = 0;

/* > \brief \b DLASY2 solves the Sylvester matrix equation where the matrices are of order 1 or 2. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASY2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasy2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasy2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasy2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, */
/*                          LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            LTRANL, LTRANR */
/*       INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2 */
/*       DOUBLE PRECISION   SCALE, XNORM */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in */
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
/* >          TL is DOUBLE PRECISION array, dimension (LDTL,2) */
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
/* >          TR is DOUBLE PRECISION array, dimension (LDTR,2) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,2) */
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
/* >          SCALE is DOUBLE PRECISION */
/* >          On exit, SCALE contains the scale factor. SCALE is chosen */
/* >          less than or equal to 1 to prevent the solution overflowing. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,2) */
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
/* >          XNORM is DOUBLE PRECISION */
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

/* > \ingroup doubleSYauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlasy2_(logical *ltranl, logical *ltranr, integer *isgn, 
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
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical xswap;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
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
#line 224 "dlasy2.f"
    /* Parameter adjustments */
#line 224 "dlasy2.f"
    tl_dim1 = *ldtl;
#line 224 "dlasy2.f"
    tl_offset = 1 + tl_dim1;
#line 224 "dlasy2.f"
    tl -= tl_offset;
#line 224 "dlasy2.f"
    tr_dim1 = *ldtr;
#line 224 "dlasy2.f"
    tr_offset = 1 + tr_dim1;
#line 224 "dlasy2.f"
    tr -= tr_offset;
#line 224 "dlasy2.f"
    b_dim1 = *ldb;
#line 224 "dlasy2.f"
    b_offset = 1 + b_dim1;
#line 224 "dlasy2.f"
    b -= b_offset;
#line 224 "dlasy2.f"
    x_dim1 = *ldx;
#line 224 "dlasy2.f"
    x_offset = 1 + x_dim1;
#line 224 "dlasy2.f"
    x -= x_offset;
#line 224 "dlasy2.f"

#line 224 "dlasy2.f"
    /* Function Body */
/*     .. */
/*     .. Executable Statements .. */

/*     Do not check the input parameters for errors */

#line 233 "dlasy2.f"
    *info = 0;

/*     Quick return if possible */

#line 237 "dlasy2.f"
    if (*n1 == 0 || *n2 == 0) {
#line 237 "dlasy2.f"
	return 0;
#line 237 "dlasy2.f"
    }

/*     Set constants to control overflow */

#line 242 "dlasy2.f"
    eps = dlamch_("P", (ftnlen)1);
#line 243 "dlasy2.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 244 "dlasy2.f"
    sgn = (doublereal) (*isgn);

#line 246 "dlasy2.f"
    k = *n1 + *n1 + *n2 - 2;
#line 247 "dlasy2.f"
    switch (k) {
#line 247 "dlasy2.f"
	case 1:  goto L10;
#line 247 "dlasy2.f"
	case 2:  goto L20;
#line 247 "dlasy2.f"
	case 3:  goto L30;
#line 247 "dlasy2.f"
	case 4:  goto L50;
#line 247 "dlasy2.f"
    }

/*     1 by 1: TL11*X + SGN*X*TR11 = B11 */

#line 251 "dlasy2.f"
L10:
#line 252 "dlasy2.f"
    tau1 = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
#line 253 "dlasy2.f"
    bet = abs(tau1);
#line 254 "dlasy2.f"
    if (bet <= smlnum) {
#line 255 "dlasy2.f"
	tau1 = smlnum;
#line 256 "dlasy2.f"
	bet = smlnum;
#line 257 "dlasy2.f"
	*info = 1;
#line 258 "dlasy2.f"
    }

#line 260 "dlasy2.f"
    *scale = 1.;
#line 261 "dlasy2.f"
    gam = (d__1 = b[b_dim1 + 1], abs(d__1));
#line 262 "dlasy2.f"
    if (smlnum * gam > bet) {
#line 262 "dlasy2.f"
	*scale = 1. / gam;
#line 262 "dlasy2.f"
    }

#line 265 "dlasy2.f"
    x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / tau1;
#line 266 "dlasy2.f"
    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1));
#line 267 "dlasy2.f"
    return 0;

/*     1 by 2: */
/*     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12] */
/*                                       [TR21 TR22] */

#line 273 "dlasy2.f"
L20:

/* Computing MAX */
/* Computing MAX */
#line 275 "dlasy2.f"
    d__7 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__8 = (d__2 = tr[tr_dim1 + 1]
	    , abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = tr[(tr_dim1 <<
	     1) + 1], abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = tr[
	    tr_dim1 + 2], abs(d__4)), d__7 = max(d__7,d__8), d__8 = (d__5 = 
	    tr[(tr_dim1 << 1) + 2], abs(d__5));
#line 275 "dlasy2.f"
    d__6 = eps * max(d__7,d__8);
#line 275 "dlasy2.f"
    smin = max(d__6,smlnum);
#line 278 "dlasy2.f"
    tmp[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
#line 279 "dlasy2.f"
    tmp[3] = tl[tl_dim1 + 1] + sgn * tr[(tr_dim1 << 1) + 2];
#line 280 "dlasy2.f"
    if (*ltranr) {
#line 281 "dlasy2.f"
	tmp[1] = sgn * tr[tr_dim1 + 2];
#line 282 "dlasy2.f"
	tmp[2] = sgn * tr[(tr_dim1 << 1) + 1];
#line 283 "dlasy2.f"
    } else {
#line 284 "dlasy2.f"
	tmp[1] = sgn * tr[(tr_dim1 << 1) + 1];
#line 285 "dlasy2.f"
	tmp[2] = sgn * tr[tr_dim1 + 2];
#line 286 "dlasy2.f"
    }
#line 287 "dlasy2.f"
    btmp[0] = b[b_dim1 + 1];
#line 288 "dlasy2.f"
    btmp[1] = b[(b_dim1 << 1) + 1];
#line 289 "dlasy2.f"
    goto L40;

/*     2 by 1: */
/*          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11] */
/*            [TL21 TL22] [X21]         [X21]         [B21] */

#line 295 "dlasy2.f"
L30:
/* Computing MAX */
/* Computing MAX */
#line 296 "dlasy2.f"
    d__7 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__8 = (d__2 = tl[tl_dim1 + 1]
	    , abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = tl[(tl_dim1 <<
	     1) + 1], abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = tl[
	    tl_dim1 + 2], abs(d__4)), d__7 = max(d__7,d__8), d__8 = (d__5 = 
	    tl[(tl_dim1 << 1) + 2], abs(d__5));
#line 296 "dlasy2.f"
    d__6 = eps * max(d__7,d__8);
#line 296 "dlasy2.f"
    smin = max(d__6,smlnum);
#line 299 "dlasy2.f"
    tmp[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
#line 300 "dlasy2.f"
    tmp[3] = tl[(tl_dim1 << 1) + 2] + sgn * tr[tr_dim1 + 1];
#line 301 "dlasy2.f"
    if (*ltranl) {
#line 302 "dlasy2.f"
	tmp[1] = tl[(tl_dim1 << 1) + 1];
#line 303 "dlasy2.f"
	tmp[2] = tl[tl_dim1 + 2];
#line 304 "dlasy2.f"
    } else {
#line 305 "dlasy2.f"
	tmp[1] = tl[tl_dim1 + 2];
#line 306 "dlasy2.f"
	tmp[2] = tl[(tl_dim1 << 1) + 1];
#line 307 "dlasy2.f"
    }
#line 308 "dlasy2.f"
    btmp[0] = b[b_dim1 + 1];
#line 309 "dlasy2.f"
    btmp[1] = b[b_dim1 + 2];
#line 310 "dlasy2.f"
L40:

/*     Solve 2 by 2 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 315 "dlasy2.f"
    ipiv = idamax_(&c__4, tmp, &c__1);
#line 316 "dlasy2.f"
    u11 = tmp[ipiv - 1];
#line 317 "dlasy2.f"
    if (abs(u11) <= smin) {
#line 318 "dlasy2.f"
	*info = 1;
#line 319 "dlasy2.f"
	u11 = smin;
#line 320 "dlasy2.f"
    }
#line 321 "dlasy2.f"
    u12 = tmp[locu12[ipiv - 1] - 1];
#line 322 "dlasy2.f"
    l21 = tmp[locl21[ipiv - 1] - 1] / u11;
#line 323 "dlasy2.f"
    u22 = tmp[locu22[ipiv - 1] - 1] - u12 * l21;
#line 324 "dlasy2.f"
    xswap = xswpiv[ipiv - 1];
#line 325 "dlasy2.f"
    bswap = bswpiv[ipiv - 1];
#line 326 "dlasy2.f"
    if (abs(u22) <= smin) {
#line 327 "dlasy2.f"
	*info = 1;
#line 328 "dlasy2.f"
	u22 = smin;
#line 329 "dlasy2.f"
    }
#line 330 "dlasy2.f"
    if (bswap) {
#line 331 "dlasy2.f"
	temp = btmp[1];
#line 332 "dlasy2.f"
	btmp[1] = btmp[0] - l21 * temp;
#line 333 "dlasy2.f"
	btmp[0] = temp;
#line 334 "dlasy2.f"
    } else {
#line 335 "dlasy2.f"
	btmp[1] -= l21 * btmp[0];
#line 336 "dlasy2.f"
    }
#line 337 "dlasy2.f"
    *scale = 1.;
#line 338 "dlasy2.f"
    if (smlnum * 2. * abs(btmp[1]) > abs(u22) || smlnum * 2. * abs(btmp[0]) > 
	    abs(u11)) {
/* Computing MAX */
#line 340 "dlasy2.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]);
#line 340 "dlasy2.f"
	*scale = .5 / max(d__1,d__2);
#line 341 "dlasy2.f"
	btmp[0] *= *scale;
#line 342 "dlasy2.f"
	btmp[1] *= *scale;
#line 343 "dlasy2.f"
    }
#line 344 "dlasy2.f"
    x2[1] = btmp[1] / u22;
#line 345 "dlasy2.f"
    x2[0] = btmp[0] / u11 - u12 / u11 * x2[1];
#line 346 "dlasy2.f"
    if (xswap) {
#line 347 "dlasy2.f"
	temp = x2[1];
#line 348 "dlasy2.f"
	x2[1] = x2[0];
#line 349 "dlasy2.f"
	x2[0] = temp;
#line 350 "dlasy2.f"
    }
#line 351 "dlasy2.f"
    x[x_dim1 + 1] = x2[0];
#line 352 "dlasy2.f"
    if (*n1 == 1) {
#line 353 "dlasy2.f"
	x[(x_dim1 << 1) + 1] = x2[1];
#line 354 "dlasy2.f"
	*xnorm = (d__1 = x[x_dim1 + 1], abs(d__1)) + (d__2 = x[(x_dim1 << 1) 
		+ 1], abs(d__2));
#line 355 "dlasy2.f"
    } else {
#line 356 "dlasy2.f"
	x[x_dim1 + 2] = x2[1];
/* Computing MAX */
#line 357 "dlasy2.f"
	d__3 = (d__1 = x[x_dim1 + 1], abs(d__1)), d__4 = (d__2 = x[x_dim1 + 2]
		, abs(d__2));
#line 357 "dlasy2.f"
	*xnorm = max(d__3,d__4);
#line 358 "dlasy2.f"
    }
#line 359 "dlasy2.f"
    return 0;

/*     2 by 2: */
/*     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12] */
/*       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22] */

/*     Solve equivalent 4 by 4 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 368 "dlasy2.f"
L50:
/* Computing MAX */
#line 369 "dlasy2.f"
    d__5 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__6 = (d__2 = tr[(tr_dim1 << 
	    1) + 1], abs(d__2)), d__5 = max(d__5,d__6), d__6 = (d__3 = tr[
	    tr_dim1 + 2], abs(d__3)), d__5 = max(d__5,d__6), d__6 = (d__4 = 
	    tr[(tr_dim1 << 1) + 2], abs(d__4));
#line 369 "dlasy2.f"
    smin = max(d__5,d__6);
/* Computing MAX */
#line 371 "dlasy2.f"
    d__5 = smin, d__6 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__5 = max(d__5,
	    d__6), d__6 = (d__2 = tl[(tl_dim1 << 1) + 1], abs(d__2)), d__5 = 
	    max(d__5,d__6), d__6 = (d__3 = tl[tl_dim1 + 2], abs(d__3)), d__5 =
	     max(d__5,d__6), d__6 = (d__4 = tl[(tl_dim1 << 1) + 2], abs(d__4))
	    ;
#line 371 "dlasy2.f"
    smin = max(d__5,d__6);
/* Computing MAX */
#line 373 "dlasy2.f"
    d__1 = eps * smin;
#line 373 "dlasy2.f"
    smin = max(d__1,smlnum);
#line 374 "dlasy2.f"
    btmp[0] = 0.;
#line 375 "dlasy2.f"
    dcopy_(&c__16, btmp, &c__0, t16, &c__1);
#line 376 "dlasy2.f"
    t16[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
#line 377 "dlasy2.f"
    t16[5] = tl[(tl_dim1 << 1) + 2] + sgn * tr[tr_dim1 + 1];
#line 378 "dlasy2.f"
    t16[10] = tl[tl_dim1 + 1] + sgn * tr[(tr_dim1 << 1) + 2];
#line 379 "dlasy2.f"
    t16[15] = tl[(tl_dim1 << 1) + 2] + sgn * tr[(tr_dim1 << 1) + 2];
#line 380 "dlasy2.f"
    if (*ltranl) {
#line 381 "dlasy2.f"
	t16[4] = tl[tl_dim1 + 2];
#line 382 "dlasy2.f"
	t16[1] = tl[(tl_dim1 << 1) + 1];
#line 383 "dlasy2.f"
	t16[14] = tl[tl_dim1 + 2];
#line 384 "dlasy2.f"
	t16[11] = tl[(tl_dim1 << 1) + 1];
#line 385 "dlasy2.f"
    } else {
#line 386 "dlasy2.f"
	t16[4] = tl[(tl_dim1 << 1) + 1];
#line 387 "dlasy2.f"
	t16[1] = tl[tl_dim1 + 2];
#line 388 "dlasy2.f"
	t16[14] = tl[(tl_dim1 << 1) + 1];
#line 389 "dlasy2.f"
	t16[11] = tl[tl_dim1 + 2];
#line 390 "dlasy2.f"
    }
#line 391 "dlasy2.f"
    if (*ltranr) {
#line 392 "dlasy2.f"
	t16[8] = sgn * tr[(tr_dim1 << 1) + 1];
#line 393 "dlasy2.f"
	t16[13] = sgn * tr[(tr_dim1 << 1) + 1];
#line 394 "dlasy2.f"
	t16[2] = sgn * tr[tr_dim1 + 2];
#line 395 "dlasy2.f"
	t16[7] = sgn * tr[tr_dim1 + 2];
#line 396 "dlasy2.f"
    } else {
#line 397 "dlasy2.f"
	t16[8] = sgn * tr[tr_dim1 + 2];
#line 398 "dlasy2.f"
	t16[13] = sgn * tr[tr_dim1 + 2];
#line 399 "dlasy2.f"
	t16[2] = sgn * tr[(tr_dim1 << 1) + 1];
#line 400 "dlasy2.f"
	t16[7] = sgn * tr[(tr_dim1 << 1) + 1];
#line 401 "dlasy2.f"
    }
#line 402 "dlasy2.f"
    btmp[0] = b[b_dim1 + 1];
#line 403 "dlasy2.f"
    btmp[1] = b[b_dim1 + 2];
#line 404 "dlasy2.f"
    btmp[2] = b[(b_dim1 << 1) + 1];
#line 405 "dlasy2.f"
    btmp[3] = b[(b_dim1 << 1) + 2];

/*     Perform elimination */

#line 409 "dlasy2.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 410 "dlasy2.f"
	xmax = 0.;
#line 411 "dlasy2.f"
	for (ip = i__; ip <= 4; ++ip) {
#line 412 "dlasy2.f"
	    for (jp = i__; jp <= 4; ++jp) {
#line 413 "dlasy2.f"
		if ((d__1 = t16[ip + (jp << 2) - 5], abs(d__1)) >= xmax) {
#line 414 "dlasy2.f"
		    xmax = (d__1 = t16[ip + (jp << 2) - 5], abs(d__1));
#line 415 "dlasy2.f"
		    ipsv = ip;
#line 416 "dlasy2.f"
		    jpsv = jp;
#line 417 "dlasy2.f"
		}
#line 418 "dlasy2.f"
/* L60: */
#line 418 "dlasy2.f"
	    }
#line 419 "dlasy2.f"
/* L70: */
#line 419 "dlasy2.f"
	}
#line 420 "dlasy2.f"
	if (ipsv != i__) {
#line 421 "dlasy2.f"
	    dswap_(&c__4, &t16[ipsv - 1], &c__4, &t16[i__ - 1], &c__4);
#line 422 "dlasy2.f"
	    temp = btmp[i__ - 1];
#line 423 "dlasy2.f"
	    btmp[i__ - 1] = btmp[ipsv - 1];
#line 424 "dlasy2.f"
	    btmp[ipsv - 1] = temp;
#line 425 "dlasy2.f"
	}
#line 426 "dlasy2.f"
	if (jpsv != i__) {
#line 426 "dlasy2.f"
	    dswap_(&c__4, &t16[(jpsv << 2) - 4], &c__1, &t16[(i__ << 2) - 4], 
		    &c__1);
#line 426 "dlasy2.f"
	}
#line 428 "dlasy2.f"
	jpiv[i__ - 1] = jpsv;
#line 429 "dlasy2.f"
	if ((d__1 = t16[i__ + (i__ << 2) - 5], abs(d__1)) < smin) {
#line 430 "dlasy2.f"
	    *info = 1;
#line 431 "dlasy2.f"
	    t16[i__ + (i__ << 2) - 5] = smin;
#line 432 "dlasy2.f"
	}
#line 433 "dlasy2.f"
	for (j = i__ + 1; j <= 4; ++j) {
#line 434 "dlasy2.f"
	    t16[j + (i__ << 2) - 5] /= t16[i__ + (i__ << 2) - 5];
#line 435 "dlasy2.f"
	    btmp[j - 1] -= t16[j + (i__ << 2) - 5] * btmp[i__ - 1];
#line 436 "dlasy2.f"
	    for (k = i__ + 1; k <= 4; ++k) {
#line 437 "dlasy2.f"
		t16[j + (k << 2) - 5] -= t16[j + (i__ << 2) - 5] * t16[i__ + (
			k << 2) - 5];
#line 438 "dlasy2.f"
/* L80: */
#line 438 "dlasy2.f"
	    }
#line 439 "dlasy2.f"
/* L90: */
#line 439 "dlasy2.f"
	}
#line 440 "dlasy2.f"
/* L100: */
#line 440 "dlasy2.f"
    }
#line 441 "dlasy2.f"
    if (abs(t16[15]) < smin) {
#line 441 "dlasy2.f"
	t16[15] = smin;
#line 441 "dlasy2.f"
    }
#line 443 "dlasy2.f"
    *scale = 1.;
#line 444 "dlasy2.f"
    if (smlnum * 8. * abs(btmp[0]) > abs(t16[0]) || smlnum * 8. * abs(btmp[1])
	     > abs(t16[5]) || smlnum * 8. * abs(btmp[2]) > abs(t16[10]) || 
	    smlnum * 8. * abs(btmp[3]) > abs(t16[15])) {
/* Computing MAX */
#line 448 "dlasy2.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]), d__1 = max(d__1,d__2), d__2 
		= abs(btmp[2]), d__1 = max(d__1,d__2), d__2 = abs(btmp[3]);
#line 448 "dlasy2.f"
	*scale = .125 / max(d__1,d__2);
#line 450 "dlasy2.f"
	btmp[0] *= *scale;
#line 451 "dlasy2.f"
	btmp[1] *= *scale;
#line 452 "dlasy2.f"
	btmp[2] *= *scale;
#line 453 "dlasy2.f"
	btmp[3] *= *scale;
#line 454 "dlasy2.f"
    }
#line 455 "dlasy2.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 456 "dlasy2.f"
	k = 5 - i__;
#line 457 "dlasy2.f"
	temp = 1. / t16[k + (k << 2) - 5];
#line 458 "dlasy2.f"
	tmp[k - 1] = btmp[k - 1] * temp;
#line 459 "dlasy2.f"
	for (j = k + 1; j <= 4; ++j) {
#line 460 "dlasy2.f"
	    tmp[k - 1] -= temp * t16[k + (j << 2) - 5] * tmp[j - 1];
#line 461 "dlasy2.f"
/* L110: */
#line 461 "dlasy2.f"
	}
#line 462 "dlasy2.f"
/* L120: */
#line 462 "dlasy2.f"
    }
#line 463 "dlasy2.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 464 "dlasy2.f"
	if (jpiv[4 - i__ - 1] != 4 - i__) {
#line 465 "dlasy2.f"
	    temp = tmp[4 - i__ - 1];
#line 466 "dlasy2.f"
	    tmp[4 - i__ - 1] = tmp[jpiv[4 - i__ - 1] - 1];
#line 467 "dlasy2.f"
	    tmp[jpiv[4 - i__ - 1] - 1] = temp;
#line 468 "dlasy2.f"
	}
#line 469 "dlasy2.f"
/* L130: */
#line 469 "dlasy2.f"
    }
#line 470 "dlasy2.f"
    x[x_dim1 + 1] = tmp[0];
#line 471 "dlasy2.f"
    x[x_dim1 + 2] = tmp[1];
#line 472 "dlasy2.f"
    x[(x_dim1 << 1) + 1] = tmp[2];
#line 473 "dlasy2.f"
    x[(x_dim1 << 1) + 2] = tmp[3];
/* Computing MAX */
#line 474 "dlasy2.f"
    d__1 = abs(tmp[0]) + abs(tmp[2]), d__2 = abs(tmp[1]) + abs(tmp[3]);
#line 474 "dlasy2.f"
    *xnorm = max(d__1,d__2);
#line 476 "dlasy2.f"
    return 0;

/*     End of DLASY2 */

} /* dlasy2_ */

