#line 1 "dlaln2.f"
/* dlaln2.f -- translated by f2c (version 20100827).
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

#line 1 "dlaln2.f"
/* > \brief \b DLALN2 solves a 1-by-1 or 2-by-2 linear system of equations of the specified form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLALN2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaln2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaln2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaln2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, */
/*                          LDB, WR, WI, X, LDX, SCALE, XNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            LTRANS */
/*       INTEGER            INFO, LDA, LDB, LDX, NA, NW */
/*       DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLALN2 solves a system of the form  (ca A - w D ) X = s B */
/* > or (ca A**T - w D) X = s B   with possible scaling ("s") and */
/* > perturbation of A.  (A**T means A-transpose.) */
/* > */
/* > A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA */
/* > real diagonal matrix, w is a real or complex value, and X and B are */
/* > NA x 1 matrices -- real if w is real, complex if w is complex.  NA */
/* > may be 1 or 2. */
/* > */
/* > If w is complex, X and B are represented as NA x 2 matrices, */
/* > the first column of each being the real part and the second */
/* > being the imaginary part. */
/* > */
/* > "s" is a scaling factor (.LE. 1), computed by DLALN2, which is */
/* > so chosen that X can be computed without overflow.  X is further */
/* > scaled if necessary to assure that norm(ca A - w D)*norm(X) is less */
/* > than overflow. */
/* > */
/* > If both singular values of (ca A - w D) are less than SMIN, */
/* > SMIN*identity will be used instead of (ca A - w D).  If only one */
/* > singular value is less than SMIN, one element of (ca A - w D) will be */
/* > perturbed enough to make the smallest singular value roughly SMIN. */
/* > If both singular values are at least SMIN, (ca A - w D) will not be */
/* > perturbed.  In any case, the perturbation will be at most some small */
/* > multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values */
/* > are computed by infinity-norm approximations, and thus will only be */
/* > correct to a factor of 2 or so. */
/* > */
/* > Note: all input quantities are assumed to be smaller than overflow */
/* > by a reasonable factor.  (See BIGNUM.) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] LTRANS */
/* > \verbatim */
/* >          LTRANS is LOGICAL */
/* >          =.TRUE.:  A-transpose will be used. */
/* >          =.FALSE.: A will be used (not transposed.) */
/* > \endverbatim */
/* > */
/* > \param[in] NA */
/* > \verbatim */
/* >          NA is INTEGER */
/* >          The size of the matrix A.  It may (only) be 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] NW */
/* > \verbatim */
/* >          NW is INTEGER */
/* >          1 if "w" is real, 2 if "w" is complex.  It may only be 1 */
/* >          or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] SMIN */
/* > \verbatim */
/* >          SMIN is DOUBLE PRECISION */
/* >          The desired lower bound on the singular values of A.  This */
/* >          should be a safe distance away from underflow or overflow, */
/* >          say, between (underflow/machine precision) and  (machine */
/* >          precision * overflow ).  (See BIGNUM and ULP.) */
/* > \endverbatim */
/* > */
/* > \param[in] CA */
/* > \verbatim */
/* >          CA is DOUBLE PRECISION */
/* >          The coefficient c, which A is multiplied by. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,NA) */
/* >          The NA x NA matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of A.  It must be at least NA. */
/* > \endverbatim */
/* > */
/* > \param[in] D1 */
/* > \verbatim */
/* >          D1 is DOUBLE PRECISION */
/* >          The 1,1 element in the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] D2 */
/* > \verbatim */
/* >          D2 is DOUBLE PRECISION */
/* >          The 2,2 element in the diagonal matrix D.  Not used if NA=1. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NW) */
/* >          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is */
/* >          complex), column 1 contains the real part of B and column 2 */
/* >          contains the imaginary part. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of B.  It must be at least NA. */
/* > \endverbatim */
/* > */
/* > \param[in] WR */
/* > \verbatim */
/* >          WR is DOUBLE PRECISION */
/* >          The real part of the scalar "w". */
/* > \endverbatim */
/* > */
/* > \param[in] WI */
/* > \verbatim */
/* >          WI is DOUBLE PRECISION */
/* >          The imaginary part of the scalar "w".  Not used if NW=1. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,NW) */
/* >          The NA x NW matrix X (unknowns), as computed by DLALN2. */
/* >          If NW=2 ("w" is complex), on exit, column 1 will contain */
/* >          the real part of X and column 2 will contain the imaginary */
/* >          part. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of X.  It must be at least NA. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          The scale factor that B must be multiplied by to insure */
/* >          that overflow does not occur when computing X.  Thus, */
/* >          (ca A - w D) X  will be SCALE*B, not B (ignoring */
/* >          perturbations of A.)  It will be at most 1. */
/* > \endverbatim */
/* > */
/* > \param[out] XNORM */
/* > \verbatim */
/* >          XNORM is DOUBLE PRECISION */
/* >          The infinity-norm of X, when X is regarded as an NA x NW */
/* >          real matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          An error flag.  It will be set to zero if no error occurs, */
/* >          a negative number if an argument is in error, or a positive */
/* >          number if  ca A - w D  had to be perturbed. */
/* >          The possible values are: */
/* >          = 0: No error occurred, and (ca A - w D) did not have to be */
/* >                 perturbed. */
/* >          = 1: (ca A - w D) had to be perturbed to make its smallest */
/* >               (or only) singular value greater than SMIN. */
/* >          NOTE: In the interests of speed, this routine does not */
/* >                check the inputs for errors. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaln2_(logical *ltrans, integer *na, integer *nw, 
	doublereal *smin, doublereal *ca, doublereal *a, integer *lda, 
	doublereal *d1, doublereal *d2, doublereal *b, integer *ldb, 
	doublereal *wr, doublereal *wi, doublereal *x, integer *ldx, 
	doublereal *scale, doublereal *xnorm, integer *info)
{
    /* Initialized data */

    static logical zswap[4] = { FALSE_,FALSE_,TRUE_,TRUE_ };
    static logical rswap[4] = { FALSE_,TRUE_,FALSE_,TRUE_ };
    static integer ipivot[16]	/* was [4][4] */ = { 1,2,3,4,2,1,4,3,3,4,1,2,
	    4,3,2,1 };

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    static doublereal equiv_0[4], equiv_1[4];

    /* Local variables */
    static integer j;
#define ci (equiv_0)
#define cr (equiv_1)
    static doublereal bi1, bi2, br1, br2, xi1, xi2, xr1, xr2, ci21, ci22, 
	    cr21, cr22, li21, csi, ui11, lr21, ui12, ui22;
#define civ (equiv_0)
    static doublereal csr, ur11, ur12, ur22;
#define crv (equiv_1)
    static doublereal bbnd, cmax, ui11r, ui12s, temp, ur11r, ur12s, u22abs;
    static integer icmax;
    static doublereal bnorm, cnorm, smini;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal bignum, smlnum;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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
/*     .. Equivalences .. */
/*     .. */
/*     .. Data statements .. */
#line 271 "dlaln2.f"
    /* Parameter adjustments */
#line 271 "dlaln2.f"
    a_dim1 = *lda;
#line 271 "dlaln2.f"
    a_offset = 1 + a_dim1;
#line 271 "dlaln2.f"
    a -= a_offset;
#line 271 "dlaln2.f"
    b_dim1 = *ldb;
#line 271 "dlaln2.f"
    b_offset = 1 + b_dim1;
#line 271 "dlaln2.f"
    b -= b_offset;
#line 271 "dlaln2.f"
    x_dim1 = *ldx;
#line 271 "dlaln2.f"
    x_offset = 1 + x_dim1;
#line 271 "dlaln2.f"
    x -= x_offset;
#line 271 "dlaln2.f"

#line 271 "dlaln2.f"
    /* Function Body */
/*     .. */
/*     .. Executable Statements .. */

/*     Compute BIGNUM */

#line 280 "dlaln2.f"
    smlnum = 2. * dlamch_("Safe minimum", (ftnlen)12);
#line 281 "dlaln2.f"
    bignum = 1. / smlnum;
#line 282 "dlaln2.f"
    smini = max(*smin,smlnum);

/*     Don't check for input errors */

#line 286 "dlaln2.f"
    *info = 0;

/*     Standard Initializations */

#line 290 "dlaln2.f"
    *scale = 1.;

#line 292 "dlaln2.f"
    if (*na == 1) {

/*        1 x 1  (i.e., scalar) system   C X = B */

#line 296 "dlaln2.f"
	if (*nw == 1) {

/*           Real 1x1 system. */

/*           C = ca A - w D */

#line 302 "dlaln2.f"
	    csr = *ca * a[a_dim1 + 1] - *wr * *d1;
#line 303 "dlaln2.f"
	    cnorm = abs(csr);

/*           If | C | < SMINI, use C = SMINI */

#line 307 "dlaln2.f"
	    if (cnorm < smini) {
#line 308 "dlaln2.f"
		csr = smini;
#line 309 "dlaln2.f"
		cnorm = smini;
#line 310 "dlaln2.f"
		*info = 1;
#line 311 "dlaln2.f"
	    }

/*           Check scaling for  X = B / C */

#line 315 "dlaln2.f"
	    bnorm = (d__1 = b[b_dim1 + 1], abs(d__1));
#line 316 "dlaln2.f"
	    if (cnorm < 1. && bnorm > 1.) {
#line 317 "dlaln2.f"
		if (bnorm > bignum * cnorm) {
#line 317 "dlaln2.f"
		    *scale = 1. / bnorm;
#line 317 "dlaln2.f"
		}
#line 319 "dlaln2.f"
	    }

/*           Compute X */

#line 323 "dlaln2.f"
	    x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / csr;
#line 324 "dlaln2.f"
	    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1));
#line 325 "dlaln2.f"
	} else {

/*           Complex 1x1 system (w is complex) */

/*           C = ca A - w D */

#line 331 "dlaln2.f"
	    csr = *ca * a[a_dim1 + 1] - *wr * *d1;
#line 332 "dlaln2.f"
	    csi = -(*wi) * *d1;
#line 333 "dlaln2.f"
	    cnorm = abs(csr) + abs(csi);

/*           If | C | < SMINI, use C = SMINI */

#line 337 "dlaln2.f"
	    if (cnorm < smini) {
#line 338 "dlaln2.f"
		csr = smini;
#line 339 "dlaln2.f"
		csi = 0.;
#line 340 "dlaln2.f"
		cnorm = smini;
#line 341 "dlaln2.f"
		*info = 1;
#line 342 "dlaln2.f"
	    }

/*           Check scaling for  X = B / C */

#line 346 "dlaln2.f"
	    bnorm = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 
		    1) + 1], abs(d__2));
#line 347 "dlaln2.f"
	    if (cnorm < 1. && bnorm > 1.) {
#line 348 "dlaln2.f"
		if (bnorm > bignum * cnorm) {
#line 348 "dlaln2.f"
		    *scale = 1. / bnorm;
#line 348 "dlaln2.f"
		}
#line 350 "dlaln2.f"
	    }

/*           Compute X */

#line 354 "dlaln2.f"
	    d__1 = *scale * b[b_dim1 + 1];
#line 354 "dlaln2.f"
	    d__2 = *scale * b[(b_dim1 << 1) + 1];
#line 354 "dlaln2.f"
	    dladiv_(&d__1, &d__2, &csr, &csi, &x[x_dim1 + 1], &x[(x_dim1 << 1)
		     + 1]);
#line 356 "dlaln2.f"
	    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1)) + (d__2 = x[(x_dim1 << 
		    1) + 1], abs(d__2));
#line 357 "dlaln2.f"
	}

#line 359 "dlaln2.f"
    } else {

/*        2x2 System */

/*        Compute the real part of  C = ca A - w D  (or  ca A**T - w D ) */

#line 365 "dlaln2.f"
	cr[0] = *ca * a[a_dim1 + 1] - *wr * *d1;
#line 366 "dlaln2.f"
	cr[3] = *ca * a[(a_dim1 << 1) + 2] - *wr * *d2;
#line 367 "dlaln2.f"
	if (*ltrans) {
#line 368 "dlaln2.f"
	    cr[2] = *ca * a[a_dim1 + 2];
#line 369 "dlaln2.f"
	    cr[1] = *ca * a[(a_dim1 << 1) + 1];
#line 370 "dlaln2.f"
	} else {
#line 371 "dlaln2.f"
	    cr[1] = *ca * a[a_dim1 + 2];
#line 372 "dlaln2.f"
	    cr[2] = *ca * a[(a_dim1 << 1) + 1];
#line 373 "dlaln2.f"
	}

#line 375 "dlaln2.f"
	if (*nw == 1) {

/*           Real 2x2 system  (w is real) */

/*           Find the largest element in C */

#line 381 "dlaln2.f"
	    cmax = 0.;
#line 382 "dlaln2.f"
	    icmax = 0;

#line 384 "dlaln2.f"
	    for (j = 1; j <= 4; ++j) {
#line 385 "dlaln2.f"
		if ((d__1 = crv[j - 1], abs(d__1)) > cmax) {
#line 386 "dlaln2.f"
		    cmax = (d__1 = crv[j - 1], abs(d__1));
#line 387 "dlaln2.f"
		    icmax = j;
#line 388 "dlaln2.f"
		}
#line 389 "dlaln2.f"
/* L10: */
#line 389 "dlaln2.f"
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

#line 393 "dlaln2.f"
	    if (cmax < smini) {
/* Computing MAX */
#line 394 "dlaln2.f"
		d__3 = (d__1 = b[b_dim1 + 1], abs(d__1)), d__4 = (d__2 = b[
			b_dim1 + 2], abs(d__2));
#line 394 "dlaln2.f"
		bnorm = max(d__3,d__4);
#line 395 "dlaln2.f"
		if (smini < 1. && bnorm > 1.) {
#line 396 "dlaln2.f"
		    if (bnorm > bignum * smini) {
#line 396 "dlaln2.f"
			*scale = 1. / bnorm;
#line 396 "dlaln2.f"
		    }
#line 398 "dlaln2.f"
		}
#line 399 "dlaln2.f"
		temp = *scale / smini;
#line 400 "dlaln2.f"
		x[x_dim1 + 1] = temp * b[b_dim1 + 1];
#line 401 "dlaln2.f"
		x[x_dim1 + 2] = temp * b[b_dim1 + 2];
#line 402 "dlaln2.f"
		*xnorm = temp * bnorm;
#line 403 "dlaln2.f"
		*info = 1;
#line 404 "dlaln2.f"
		return 0;
#line 405 "dlaln2.f"
	    }

/*           Gaussian elimination with complete pivoting. */

#line 409 "dlaln2.f"
	    ur11 = crv[icmax - 1];
#line 410 "dlaln2.f"
	    cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
#line 411 "dlaln2.f"
	    ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
#line 412 "dlaln2.f"
	    cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
#line 413 "dlaln2.f"
	    ur11r = 1. / ur11;
#line 414 "dlaln2.f"
	    lr21 = ur11r * cr21;
#line 415 "dlaln2.f"
	    ur22 = cr22 - ur12 * lr21;

/*           If smaller pivot < SMINI, use SMINI */

#line 419 "dlaln2.f"
	    if (abs(ur22) < smini) {
#line 420 "dlaln2.f"
		ur22 = smini;
#line 421 "dlaln2.f"
		*info = 1;
#line 422 "dlaln2.f"
	    }
#line 423 "dlaln2.f"
	    if (rswap[icmax - 1]) {
#line 424 "dlaln2.f"
		br1 = b[b_dim1 + 2];
#line 425 "dlaln2.f"
		br2 = b[b_dim1 + 1];
#line 426 "dlaln2.f"
	    } else {
#line 427 "dlaln2.f"
		br1 = b[b_dim1 + 1];
#line 428 "dlaln2.f"
		br2 = b[b_dim1 + 2];
#line 429 "dlaln2.f"
	    }
#line 430 "dlaln2.f"
	    br2 -= lr21 * br1;
/* Computing MAX */
#line 431 "dlaln2.f"
	    d__2 = (d__1 = br1 * (ur22 * ur11r), abs(d__1)), d__3 = abs(br2);
#line 431 "dlaln2.f"
	    bbnd = max(d__2,d__3);
#line 432 "dlaln2.f"
	    if (bbnd > 1. && abs(ur22) < 1.) {
#line 433 "dlaln2.f"
		if (bbnd >= bignum * abs(ur22)) {
#line 433 "dlaln2.f"
		    *scale = 1. / bbnd;
#line 433 "dlaln2.f"
		}
#line 435 "dlaln2.f"
	    }

#line 437 "dlaln2.f"
	    xr2 = br2 * *scale / ur22;
#line 438 "dlaln2.f"
	    xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
#line 439 "dlaln2.f"
	    if (zswap[icmax - 1]) {
#line 440 "dlaln2.f"
		x[x_dim1 + 1] = xr2;
#line 441 "dlaln2.f"
		x[x_dim1 + 2] = xr1;
#line 442 "dlaln2.f"
	    } else {
#line 443 "dlaln2.f"
		x[x_dim1 + 1] = xr1;
#line 444 "dlaln2.f"
		x[x_dim1 + 2] = xr2;
#line 445 "dlaln2.f"
	    }
/* Computing MAX */
#line 446 "dlaln2.f"
	    d__1 = abs(xr1), d__2 = abs(xr2);
#line 446 "dlaln2.f"
	    *xnorm = max(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

#line 450 "dlaln2.f"
	    if (*xnorm > 1. && cmax > 1.) {
#line 451 "dlaln2.f"
		if (*xnorm > bignum / cmax) {
#line 452 "dlaln2.f"
		    temp = cmax / bignum;
#line 453 "dlaln2.f"
		    x[x_dim1 + 1] = temp * x[x_dim1 + 1];
#line 454 "dlaln2.f"
		    x[x_dim1 + 2] = temp * x[x_dim1 + 2];
#line 455 "dlaln2.f"
		    *xnorm = temp * *xnorm;
#line 456 "dlaln2.f"
		    *scale = temp * *scale;
#line 457 "dlaln2.f"
		}
#line 458 "dlaln2.f"
	    }
#line 459 "dlaln2.f"
	} else {

/*           Complex 2x2 system  (w is complex) */

/*           Find the largest element in C */

#line 465 "dlaln2.f"
	    ci[0] = -(*wi) * *d1;
#line 466 "dlaln2.f"
	    ci[1] = 0.;
#line 467 "dlaln2.f"
	    ci[2] = 0.;
#line 468 "dlaln2.f"
	    ci[3] = -(*wi) * *d2;
#line 469 "dlaln2.f"
	    cmax = 0.;
#line 470 "dlaln2.f"
	    icmax = 0;

#line 472 "dlaln2.f"
	    for (j = 1; j <= 4; ++j) {
#line 473 "dlaln2.f"
		if ((d__1 = crv[j - 1], abs(d__1)) + (d__2 = civ[j - 1], abs(
			d__2)) > cmax) {
#line 474 "dlaln2.f"
		    cmax = (d__1 = crv[j - 1], abs(d__1)) + (d__2 = civ[j - 1]
			    , abs(d__2));
#line 475 "dlaln2.f"
		    icmax = j;
#line 476 "dlaln2.f"
		}
#line 477 "dlaln2.f"
/* L20: */
#line 477 "dlaln2.f"
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

#line 481 "dlaln2.f"
	    if (cmax < smini) {
/* Computing MAX */
#line 482 "dlaln2.f"
		d__5 = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 
			<< 1) + 1], abs(d__2)), d__6 = (d__3 = b[b_dim1 + 2], 
			abs(d__3)) + (d__4 = b[(b_dim1 << 1) + 2], abs(d__4));
#line 482 "dlaln2.f"
		bnorm = max(d__5,d__6);
#line 484 "dlaln2.f"
		if (smini < 1. && bnorm > 1.) {
#line 485 "dlaln2.f"
		    if (bnorm > bignum * smini) {
#line 485 "dlaln2.f"
			*scale = 1. / bnorm;
#line 485 "dlaln2.f"
		    }
#line 487 "dlaln2.f"
		}
#line 488 "dlaln2.f"
		temp = *scale / smini;
#line 489 "dlaln2.f"
		x[x_dim1 + 1] = temp * b[b_dim1 + 1];
#line 490 "dlaln2.f"
		x[x_dim1 + 2] = temp * b[b_dim1 + 2];
#line 491 "dlaln2.f"
		x[(x_dim1 << 1) + 1] = temp * b[(b_dim1 << 1) + 1];
#line 492 "dlaln2.f"
		x[(x_dim1 << 1) + 2] = temp * b[(b_dim1 << 1) + 2];
#line 493 "dlaln2.f"
		*xnorm = temp * bnorm;
#line 494 "dlaln2.f"
		*info = 1;
#line 495 "dlaln2.f"
		return 0;
#line 496 "dlaln2.f"
	    }

/*           Gaussian elimination with complete pivoting. */

#line 500 "dlaln2.f"
	    ur11 = crv[icmax - 1];
#line 501 "dlaln2.f"
	    ui11 = civ[icmax - 1];
#line 502 "dlaln2.f"
	    cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
#line 503 "dlaln2.f"
	    ci21 = civ[ipivot[(icmax << 2) - 3] - 1];
#line 504 "dlaln2.f"
	    ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
#line 505 "dlaln2.f"
	    ui12 = civ[ipivot[(icmax << 2) - 2] - 1];
#line 506 "dlaln2.f"
	    cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
#line 507 "dlaln2.f"
	    ci22 = civ[ipivot[(icmax << 2) - 1] - 1];
#line 508 "dlaln2.f"
	    if (icmax == 1 || icmax == 4) {

/*              Code when off-diagonals of pivoted C are real */

#line 512 "dlaln2.f"
		if (abs(ur11) > abs(ui11)) {
#line 513 "dlaln2.f"
		    temp = ui11 / ur11;
/* Computing 2nd power */
#line 514 "dlaln2.f"
		    d__1 = temp;
#line 514 "dlaln2.f"
		    ur11r = 1. / (ur11 * (d__1 * d__1 + 1.));
#line 515 "dlaln2.f"
		    ui11r = -temp * ur11r;
#line 516 "dlaln2.f"
		} else {
#line 517 "dlaln2.f"
		    temp = ur11 / ui11;
/* Computing 2nd power */
#line 518 "dlaln2.f"
		    d__1 = temp;
#line 518 "dlaln2.f"
		    ui11r = -1. / (ui11 * (d__1 * d__1 + 1.));
#line 519 "dlaln2.f"
		    ur11r = -temp * ui11r;
#line 520 "dlaln2.f"
		}
#line 521 "dlaln2.f"
		lr21 = cr21 * ur11r;
#line 522 "dlaln2.f"
		li21 = cr21 * ui11r;
#line 523 "dlaln2.f"
		ur12s = ur12 * ur11r;
#line 524 "dlaln2.f"
		ui12s = ur12 * ui11r;
#line 525 "dlaln2.f"
		ur22 = cr22 - ur12 * lr21;
#line 526 "dlaln2.f"
		ui22 = ci22 - ur12 * li21;
#line 527 "dlaln2.f"
	    } else {

/*              Code when diagonals of pivoted C are real */

#line 531 "dlaln2.f"
		ur11r = 1. / ur11;
#line 532 "dlaln2.f"
		ui11r = 0.;
#line 533 "dlaln2.f"
		lr21 = cr21 * ur11r;
#line 534 "dlaln2.f"
		li21 = ci21 * ur11r;
#line 535 "dlaln2.f"
		ur12s = ur12 * ur11r;
#line 536 "dlaln2.f"
		ui12s = ui12 * ur11r;
#line 537 "dlaln2.f"
		ur22 = cr22 - ur12 * lr21 + ui12 * li21;
#line 538 "dlaln2.f"
		ui22 = -ur12 * li21 - ui12 * lr21;
#line 539 "dlaln2.f"
	    }
#line 540 "dlaln2.f"
	    u22abs = abs(ur22) + abs(ui22);

/*           If smaller pivot < SMINI, use SMINI */

#line 544 "dlaln2.f"
	    if (u22abs < smini) {
#line 545 "dlaln2.f"
		ur22 = smini;
#line 546 "dlaln2.f"
		ui22 = 0.;
#line 547 "dlaln2.f"
		*info = 1;
#line 548 "dlaln2.f"
	    }
#line 549 "dlaln2.f"
	    if (rswap[icmax - 1]) {
#line 550 "dlaln2.f"
		br2 = b[b_dim1 + 1];
#line 551 "dlaln2.f"
		br1 = b[b_dim1 + 2];
#line 552 "dlaln2.f"
		bi2 = b[(b_dim1 << 1) + 1];
#line 553 "dlaln2.f"
		bi1 = b[(b_dim1 << 1) + 2];
#line 554 "dlaln2.f"
	    } else {
#line 555 "dlaln2.f"
		br1 = b[b_dim1 + 1];
#line 556 "dlaln2.f"
		br2 = b[b_dim1 + 2];
#line 557 "dlaln2.f"
		bi1 = b[(b_dim1 << 1) + 1];
#line 558 "dlaln2.f"
		bi2 = b[(b_dim1 << 1) + 2];
#line 559 "dlaln2.f"
	    }
#line 560 "dlaln2.f"
	    br2 = br2 - lr21 * br1 + li21 * bi1;
#line 561 "dlaln2.f"
	    bi2 = bi2 - li21 * br1 - lr21 * bi1;
/* Computing MAX */
#line 562 "dlaln2.f"
	    d__1 = (abs(br1) + abs(bi1)) * (u22abs * (abs(ur11r) + abs(ui11r))
		    ), d__2 = abs(br2) + abs(bi2);
#line 562 "dlaln2.f"
	    bbnd = max(d__1,d__2);
#line 565 "dlaln2.f"
	    if (bbnd > 1. && u22abs < 1.) {
#line 566 "dlaln2.f"
		if (bbnd >= bignum * u22abs) {
#line 567 "dlaln2.f"
		    *scale = 1. / bbnd;
#line 568 "dlaln2.f"
		    br1 = *scale * br1;
#line 569 "dlaln2.f"
		    bi1 = *scale * bi1;
#line 570 "dlaln2.f"
		    br2 = *scale * br2;
#line 571 "dlaln2.f"
		    bi2 = *scale * bi2;
#line 572 "dlaln2.f"
		}
#line 573 "dlaln2.f"
	    }

#line 575 "dlaln2.f"
	    dladiv_(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
#line 576 "dlaln2.f"
	    xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
#line 577 "dlaln2.f"
	    xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
#line 578 "dlaln2.f"
	    if (zswap[icmax - 1]) {
#line 579 "dlaln2.f"
		x[x_dim1 + 1] = xr2;
#line 580 "dlaln2.f"
		x[x_dim1 + 2] = xr1;
#line 581 "dlaln2.f"
		x[(x_dim1 << 1) + 1] = xi2;
#line 582 "dlaln2.f"
		x[(x_dim1 << 1) + 2] = xi1;
#line 583 "dlaln2.f"
	    } else {
#line 584 "dlaln2.f"
		x[x_dim1 + 1] = xr1;
#line 585 "dlaln2.f"
		x[x_dim1 + 2] = xr2;
#line 586 "dlaln2.f"
		x[(x_dim1 << 1) + 1] = xi1;
#line 587 "dlaln2.f"
		x[(x_dim1 << 1) + 2] = xi2;
#line 588 "dlaln2.f"
	    }
/* Computing MAX */
#line 589 "dlaln2.f"
	    d__1 = abs(xr1) + abs(xi1), d__2 = abs(xr2) + abs(xi2);
#line 589 "dlaln2.f"
	    *xnorm = max(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

#line 593 "dlaln2.f"
	    if (*xnorm > 1. && cmax > 1.) {
#line 594 "dlaln2.f"
		if (*xnorm > bignum / cmax) {
#line 595 "dlaln2.f"
		    temp = cmax / bignum;
#line 596 "dlaln2.f"
		    x[x_dim1 + 1] = temp * x[x_dim1 + 1];
#line 597 "dlaln2.f"
		    x[x_dim1 + 2] = temp * x[x_dim1 + 2];
#line 598 "dlaln2.f"
		    x[(x_dim1 << 1) + 1] = temp * x[(x_dim1 << 1) + 1];
#line 599 "dlaln2.f"
		    x[(x_dim1 << 1) + 2] = temp * x[(x_dim1 << 1) + 2];
#line 600 "dlaln2.f"
		    *xnorm = temp * *xnorm;
#line 601 "dlaln2.f"
		    *scale = temp * *scale;
#line 602 "dlaln2.f"
		}
#line 603 "dlaln2.f"
	    }
#line 604 "dlaln2.f"
	}
#line 605 "dlaln2.f"
    }

#line 607 "dlaln2.f"
    return 0;

/*     End of DLALN2 */

} /* dlaln2_ */

#undef crv
#undef civ
#undef cr
#undef ci


