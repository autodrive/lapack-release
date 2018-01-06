#line 1 "slaln2.f"
/* slaln2.f -- translated by f2c (version 20100827).
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

#line 1 "slaln2.f"
/* > \brief \b SLALN2 solves a 1-by-1 or 2-by-2 linear system of equations of the specified form. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLALN2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaln2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaln2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaln2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, */
/*                          LDB, WR, WI, X, LDX, SCALE, XNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            LTRANS */
/*       INTEGER            INFO, LDA, LDB, LDX, NA, NW */
/*       REAL               CA, D1, D2, SCALE, SMIN, WI, WR, XNORM */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLALN2 solves a system of the form  (ca A - w D ) X = s B */
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
/* > "s" is a scaling factor (.LE. 1), computed by SLALN2, which is */
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
/* >          SMIN is REAL */
/* >          The desired lower bound on the singular values of A.  This */
/* >          should be a safe distance away from underflow or overflow, */
/* >          say, between (underflow/machine precision) and  (machine */
/* >          precision * overflow ).  (See BIGNUM and ULP.) */
/* > \endverbatim */
/* > */
/* > \param[in] CA */
/* > \verbatim */
/* >          CA is REAL */
/* >          The coefficient c, which A is multiplied by. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,NA) */
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
/* >          D1 is REAL */
/* >          The 1,1 element in the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] D2 */
/* > \verbatim */
/* >          D2 is REAL */
/* >          The 2,2 element in the diagonal matrix D.  Not used if NA=1. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NW) */
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
/* >          WR is REAL */
/* >          The real part of the scalar "w". */
/* > \endverbatim */
/* > */
/* > \param[in] WI */
/* > \verbatim */
/* >          WI is REAL */
/* >          The imaginary part of the scalar "w".  Not used if NW=1. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (LDX,NW) */
/* >          The NA x NW matrix X (unknowns), as computed by SLALN2. */
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
/* >          SCALE is REAL */
/* >          The scale factor that B must be multiplied by to insure */
/* >          that overflow does not occur when computing X.  Thus, */
/* >          (ca A - w D) X  will be SCALE*B, not B (ignoring */
/* >          perturbations of A.)  It will be at most 1. */
/* > \endverbatim */
/* > */
/* > \param[out] XNORM */
/* > \verbatim */
/* >          XNORM is REAL */
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

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slaln2_(logical *ltrans, integer *na, integer *nw, 
	doublereal *smin, doublereal *ca, doublereal *a, integer *lda, 
	doublereal *d1, doublereal *d2, doublereal *b, integer *ldb, 
	doublereal *wr, doublereal *wi, doublereal *x, integer *ldx, 
	doublereal *scale, doublereal *xnorm, integer *info)
{
    /* Initialized data */

    static logical cswap[4] = { FALSE_,FALSE_,TRUE_,TRUE_ };
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
    extern doublereal slamch_(char *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int sladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal smlnum;


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
#line 271 "slaln2.f"
    /* Parameter adjustments */
#line 271 "slaln2.f"
    a_dim1 = *lda;
#line 271 "slaln2.f"
    a_offset = 1 + a_dim1;
#line 271 "slaln2.f"
    a -= a_offset;
#line 271 "slaln2.f"
    b_dim1 = *ldb;
#line 271 "slaln2.f"
    b_offset = 1 + b_dim1;
#line 271 "slaln2.f"
    b -= b_offset;
#line 271 "slaln2.f"
    x_dim1 = *ldx;
#line 271 "slaln2.f"
    x_offset = 1 + x_dim1;
#line 271 "slaln2.f"
    x -= x_offset;
#line 271 "slaln2.f"

#line 271 "slaln2.f"
    /* Function Body */
/*     .. */
/*     .. Executable Statements .. */

/*     Compute BIGNUM */

#line 280 "slaln2.f"
    smlnum = 2. * slamch_("Safe minimum", (ftnlen)12);
#line 281 "slaln2.f"
    bignum = 1. / smlnum;
#line 282 "slaln2.f"
    smini = max(*smin,smlnum);

/*     Don't check for input errors */

#line 286 "slaln2.f"
    *info = 0;

/*     Standard Initializations */

#line 290 "slaln2.f"
    *scale = 1.;

#line 292 "slaln2.f"
    if (*na == 1) {

/*        1 x 1  (i.e., scalar) system   C X = B */

#line 296 "slaln2.f"
	if (*nw == 1) {

/*           Real 1x1 system. */

/*           C = ca A - w D */

#line 302 "slaln2.f"
	    csr = *ca * a[a_dim1 + 1] - *wr * *d1;
#line 303 "slaln2.f"
	    cnorm = abs(csr);

/*           If | C | < SMINI, use C = SMINI */

#line 307 "slaln2.f"
	    if (cnorm < smini) {
#line 308 "slaln2.f"
		csr = smini;
#line 309 "slaln2.f"
		cnorm = smini;
#line 310 "slaln2.f"
		*info = 1;
#line 311 "slaln2.f"
	    }

/*           Check scaling for  X = B / C */

#line 315 "slaln2.f"
	    bnorm = (d__1 = b[b_dim1 + 1], abs(d__1));
#line 316 "slaln2.f"
	    if (cnorm < 1. && bnorm > 1.) {
#line 317 "slaln2.f"
		if (bnorm > bignum * cnorm) {
#line 317 "slaln2.f"
		    *scale = 1. / bnorm;
#line 317 "slaln2.f"
		}
#line 319 "slaln2.f"
	    }

/*           Compute X */

#line 323 "slaln2.f"
	    x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / csr;
#line 324 "slaln2.f"
	    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1));
#line 325 "slaln2.f"
	} else {

/*           Complex 1x1 system (w is complex) */

/*           C = ca A - w D */

#line 331 "slaln2.f"
	    csr = *ca * a[a_dim1 + 1] - *wr * *d1;
#line 332 "slaln2.f"
	    csi = -(*wi) * *d1;
#line 333 "slaln2.f"
	    cnorm = abs(csr) + abs(csi);

/*           If | C | < SMINI, use C = SMINI */

#line 337 "slaln2.f"
	    if (cnorm < smini) {
#line 338 "slaln2.f"
		csr = smini;
#line 339 "slaln2.f"
		csi = 0.;
#line 340 "slaln2.f"
		cnorm = smini;
#line 341 "slaln2.f"
		*info = 1;
#line 342 "slaln2.f"
	    }

/*           Check scaling for  X = B / C */

#line 346 "slaln2.f"
	    bnorm = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 
		    1) + 1], abs(d__2));
#line 347 "slaln2.f"
	    if (cnorm < 1. && bnorm > 1.) {
#line 348 "slaln2.f"
		if (bnorm > bignum * cnorm) {
#line 348 "slaln2.f"
		    *scale = 1. / bnorm;
#line 348 "slaln2.f"
		}
#line 350 "slaln2.f"
	    }

/*           Compute X */

#line 354 "slaln2.f"
	    d__1 = *scale * b[b_dim1 + 1];
#line 354 "slaln2.f"
	    d__2 = *scale * b[(b_dim1 << 1) + 1];
#line 354 "slaln2.f"
	    sladiv_(&d__1, &d__2, &csr, &csi, &x[x_dim1 + 1], &x[(x_dim1 << 1)
		     + 1]);
#line 356 "slaln2.f"
	    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1)) + (d__2 = x[(x_dim1 << 
		    1) + 1], abs(d__2));
#line 357 "slaln2.f"
	}

#line 359 "slaln2.f"
    } else {

/*        2x2 System */

/*        Compute the real part of  C = ca A - w D  (or  ca A**T - w D ) */

#line 365 "slaln2.f"
	cr[0] = *ca * a[a_dim1 + 1] - *wr * *d1;
#line 366 "slaln2.f"
	cr[3] = *ca * a[(a_dim1 << 1) + 2] - *wr * *d2;
#line 367 "slaln2.f"
	if (*ltrans) {
#line 368 "slaln2.f"
	    cr[2] = *ca * a[a_dim1 + 2];
#line 369 "slaln2.f"
	    cr[1] = *ca * a[(a_dim1 << 1) + 1];
#line 370 "slaln2.f"
	} else {
#line 371 "slaln2.f"
	    cr[1] = *ca * a[a_dim1 + 2];
#line 372 "slaln2.f"
	    cr[2] = *ca * a[(a_dim1 << 1) + 1];
#line 373 "slaln2.f"
	}

#line 375 "slaln2.f"
	if (*nw == 1) {

/*           Real 2x2 system  (w is real) */

/*           Find the largest element in C */

#line 381 "slaln2.f"
	    cmax = 0.;
#line 382 "slaln2.f"
	    icmax = 0;

#line 384 "slaln2.f"
	    for (j = 1; j <= 4; ++j) {
#line 385 "slaln2.f"
		if ((d__1 = crv[j - 1], abs(d__1)) > cmax) {
#line 386 "slaln2.f"
		    cmax = (d__1 = crv[j - 1], abs(d__1));
#line 387 "slaln2.f"
		    icmax = j;
#line 388 "slaln2.f"
		}
#line 389 "slaln2.f"
/* L10: */
#line 389 "slaln2.f"
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

#line 393 "slaln2.f"
	    if (cmax < smini) {
/* Computing MAX */
#line 394 "slaln2.f"
		d__3 = (d__1 = b[b_dim1 + 1], abs(d__1)), d__4 = (d__2 = b[
			b_dim1 + 2], abs(d__2));
#line 394 "slaln2.f"
		bnorm = max(d__3,d__4);
#line 395 "slaln2.f"
		if (smini < 1. && bnorm > 1.) {
#line 396 "slaln2.f"
		    if (bnorm > bignum * smini) {
#line 396 "slaln2.f"
			*scale = 1. / bnorm;
#line 396 "slaln2.f"
		    }
#line 398 "slaln2.f"
		}
#line 399 "slaln2.f"
		temp = *scale / smini;
#line 400 "slaln2.f"
		x[x_dim1 + 1] = temp * b[b_dim1 + 1];
#line 401 "slaln2.f"
		x[x_dim1 + 2] = temp * b[b_dim1 + 2];
#line 402 "slaln2.f"
		*xnorm = temp * bnorm;
#line 403 "slaln2.f"
		*info = 1;
#line 404 "slaln2.f"
		return 0;
#line 405 "slaln2.f"
	    }

/*           Gaussian elimination with complete pivoting. */

#line 409 "slaln2.f"
	    ur11 = crv[icmax - 1];
#line 410 "slaln2.f"
	    cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
#line 411 "slaln2.f"
	    ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
#line 412 "slaln2.f"
	    cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
#line 413 "slaln2.f"
	    ur11r = 1. / ur11;
#line 414 "slaln2.f"
	    lr21 = ur11r * cr21;
#line 415 "slaln2.f"
	    ur22 = cr22 - ur12 * lr21;

/*           If smaller pivot < SMINI, use SMINI */

#line 419 "slaln2.f"
	    if (abs(ur22) < smini) {
#line 420 "slaln2.f"
		ur22 = smini;
#line 421 "slaln2.f"
		*info = 1;
#line 422 "slaln2.f"
	    }
#line 423 "slaln2.f"
	    if (rswap[icmax - 1]) {
#line 424 "slaln2.f"
		br1 = b[b_dim1 + 2];
#line 425 "slaln2.f"
		br2 = b[b_dim1 + 1];
#line 426 "slaln2.f"
	    } else {
#line 427 "slaln2.f"
		br1 = b[b_dim1 + 1];
#line 428 "slaln2.f"
		br2 = b[b_dim1 + 2];
#line 429 "slaln2.f"
	    }
#line 430 "slaln2.f"
	    br2 -= lr21 * br1;
/* Computing MAX */
#line 431 "slaln2.f"
	    d__2 = (d__1 = br1 * (ur22 * ur11r), abs(d__1)), d__3 = abs(br2);
#line 431 "slaln2.f"
	    bbnd = max(d__2,d__3);
#line 432 "slaln2.f"
	    if (bbnd > 1. && abs(ur22) < 1.) {
#line 433 "slaln2.f"
		if (bbnd >= bignum * abs(ur22)) {
#line 433 "slaln2.f"
		    *scale = 1. / bbnd;
#line 433 "slaln2.f"
		}
#line 435 "slaln2.f"
	    }

#line 437 "slaln2.f"
	    xr2 = br2 * *scale / ur22;
#line 438 "slaln2.f"
	    xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
#line 439 "slaln2.f"
	    if (cswap[icmax - 1]) {
#line 440 "slaln2.f"
		x[x_dim1 + 1] = xr2;
#line 441 "slaln2.f"
		x[x_dim1 + 2] = xr1;
#line 442 "slaln2.f"
	    } else {
#line 443 "slaln2.f"
		x[x_dim1 + 1] = xr1;
#line 444 "slaln2.f"
		x[x_dim1 + 2] = xr2;
#line 445 "slaln2.f"
	    }
/* Computing MAX */
#line 446 "slaln2.f"
	    d__1 = abs(xr1), d__2 = abs(xr2);
#line 446 "slaln2.f"
	    *xnorm = max(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

#line 450 "slaln2.f"
	    if (*xnorm > 1. && cmax > 1.) {
#line 451 "slaln2.f"
		if (*xnorm > bignum / cmax) {
#line 452 "slaln2.f"
		    temp = cmax / bignum;
#line 453 "slaln2.f"
		    x[x_dim1 + 1] = temp * x[x_dim1 + 1];
#line 454 "slaln2.f"
		    x[x_dim1 + 2] = temp * x[x_dim1 + 2];
#line 455 "slaln2.f"
		    *xnorm = temp * *xnorm;
#line 456 "slaln2.f"
		    *scale = temp * *scale;
#line 457 "slaln2.f"
		}
#line 458 "slaln2.f"
	    }
#line 459 "slaln2.f"
	} else {

/*           Complex 2x2 system  (w is complex) */

/*           Find the largest element in C */

#line 465 "slaln2.f"
	    ci[0] = -(*wi) * *d1;
#line 466 "slaln2.f"
	    ci[1] = 0.;
#line 467 "slaln2.f"
	    ci[2] = 0.;
#line 468 "slaln2.f"
	    ci[3] = -(*wi) * *d2;
#line 469 "slaln2.f"
	    cmax = 0.;
#line 470 "slaln2.f"
	    icmax = 0;

#line 472 "slaln2.f"
	    for (j = 1; j <= 4; ++j) {
#line 473 "slaln2.f"
		if ((d__1 = crv[j - 1], abs(d__1)) + (d__2 = civ[j - 1], abs(
			d__2)) > cmax) {
#line 474 "slaln2.f"
		    cmax = (d__1 = crv[j - 1], abs(d__1)) + (d__2 = civ[j - 1]
			    , abs(d__2));
#line 475 "slaln2.f"
		    icmax = j;
#line 476 "slaln2.f"
		}
#line 477 "slaln2.f"
/* L20: */
#line 477 "slaln2.f"
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

#line 481 "slaln2.f"
	    if (cmax < smini) {
/* Computing MAX */
#line 482 "slaln2.f"
		d__5 = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 
			<< 1) + 1], abs(d__2)), d__6 = (d__3 = b[b_dim1 + 2], 
			abs(d__3)) + (d__4 = b[(b_dim1 << 1) + 2], abs(d__4));
#line 482 "slaln2.f"
		bnorm = max(d__5,d__6);
#line 484 "slaln2.f"
		if (smini < 1. && bnorm > 1.) {
#line 485 "slaln2.f"
		    if (bnorm > bignum * smini) {
#line 485 "slaln2.f"
			*scale = 1. / bnorm;
#line 485 "slaln2.f"
		    }
#line 487 "slaln2.f"
		}
#line 488 "slaln2.f"
		temp = *scale / smini;
#line 489 "slaln2.f"
		x[x_dim1 + 1] = temp * b[b_dim1 + 1];
#line 490 "slaln2.f"
		x[x_dim1 + 2] = temp * b[b_dim1 + 2];
#line 491 "slaln2.f"
		x[(x_dim1 << 1) + 1] = temp * b[(b_dim1 << 1) + 1];
#line 492 "slaln2.f"
		x[(x_dim1 << 1) + 2] = temp * b[(b_dim1 << 1) + 2];
#line 493 "slaln2.f"
		*xnorm = temp * bnorm;
#line 494 "slaln2.f"
		*info = 1;
#line 495 "slaln2.f"
		return 0;
#line 496 "slaln2.f"
	    }

/*           Gaussian elimination with complete pivoting. */

#line 500 "slaln2.f"
	    ur11 = crv[icmax - 1];
#line 501 "slaln2.f"
	    ui11 = civ[icmax - 1];
#line 502 "slaln2.f"
	    cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
#line 503 "slaln2.f"
	    ci21 = civ[ipivot[(icmax << 2) - 3] - 1];
#line 504 "slaln2.f"
	    ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
#line 505 "slaln2.f"
	    ui12 = civ[ipivot[(icmax << 2) - 2] - 1];
#line 506 "slaln2.f"
	    cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
#line 507 "slaln2.f"
	    ci22 = civ[ipivot[(icmax << 2) - 1] - 1];
#line 508 "slaln2.f"
	    if (icmax == 1 || icmax == 4) {

/*              Code when off-diagonals of pivoted C are real */

#line 512 "slaln2.f"
		if (abs(ur11) > abs(ui11)) {
#line 513 "slaln2.f"
		    temp = ui11 / ur11;
/* Computing 2nd power */
#line 514 "slaln2.f"
		    d__1 = temp;
#line 514 "slaln2.f"
		    ur11r = 1. / (ur11 * (d__1 * d__1 + 1.));
#line 515 "slaln2.f"
		    ui11r = -temp * ur11r;
#line 516 "slaln2.f"
		} else {
#line 517 "slaln2.f"
		    temp = ur11 / ui11;
/* Computing 2nd power */
#line 518 "slaln2.f"
		    d__1 = temp;
#line 518 "slaln2.f"
		    ui11r = -1. / (ui11 * (d__1 * d__1 + 1.));
#line 519 "slaln2.f"
		    ur11r = -temp * ui11r;
#line 520 "slaln2.f"
		}
#line 521 "slaln2.f"
		lr21 = cr21 * ur11r;
#line 522 "slaln2.f"
		li21 = cr21 * ui11r;
#line 523 "slaln2.f"
		ur12s = ur12 * ur11r;
#line 524 "slaln2.f"
		ui12s = ur12 * ui11r;
#line 525 "slaln2.f"
		ur22 = cr22 - ur12 * lr21;
#line 526 "slaln2.f"
		ui22 = ci22 - ur12 * li21;
#line 527 "slaln2.f"
	    } else {

/*              Code when diagonals of pivoted C are real */

#line 531 "slaln2.f"
		ur11r = 1. / ur11;
#line 532 "slaln2.f"
		ui11r = 0.;
#line 533 "slaln2.f"
		lr21 = cr21 * ur11r;
#line 534 "slaln2.f"
		li21 = ci21 * ur11r;
#line 535 "slaln2.f"
		ur12s = ur12 * ur11r;
#line 536 "slaln2.f"
		ui12s = ui12 * ur11r;
#line 537 "slaln2.f"
		ur22 = cr22 - ur12 * lr21 + ui12 * li21;
#line 538 "slaln2.f"
		ui22 = -ur12 * li21 - ui12 * lr21;
#line 539 "slaln2.f"
	    }
#line 540 "slaln2.f"
	    u22abs = abs(ur22) + abs(ui22);

/*           If smaller pivot < SMINI, use SMINI */

#line 544 "slaln2.f"
	    if (u22abs < smini) {
#line 545 "slaln2.f"
		ur22 = smini;
#line 546 "slaln2.f"
		ui22 = 0.;
#line 547 "slaln2.f"
		*info = 1;
#line 548 "slaln2.f"
	    }
#line 549 "slaln2.f"
	    if (rswap[icmax - 1]) {
#line 550 "slaln2.f"
		br2 = b[b_dim1 + 1];
#line 551 "slaln2.f"
		br1 = b[b_dim1 + 2];
#line 552 "slaln2.f"
		bi2 = b[(b_dim1 << 1) + 1];
#line 553 "slaln2.f"
		bi1 = b[(b_dim1 << 1) + 2];
#line 554 "slaln2.f"
	    } else {
#line 555 "slaln2.f"
		br1 = b[b_dim1 + 1];
#line 556 "slaln2.f"
		br2 = b[b_dim1 + 2];
#line 557 "slaln2.f"
		bi1 = b[(b_dim1 << 1) + 1];
#line 558 "slaln2.f"
		bi2 = b[(b_dim1 << 1) + 2];
#line 559 "slaln2.f"
	    }
#line 560 "slaln2.f"
	    br2 = br2 - lr21 * br1 + li21 * bi1;
#line 561 "slaln2.f"
	    bi2 = bi2 - li21 * br1 - lr21 * bi1;
/* Computing MAX */
#line 562 "slaln2.f"
	    d__1 = (abs(br1) + abs(bi1)) * (u22abs * (abs(ur11r) + abs(ui11r))
		    ), d__2 = abs(br2) + abs(bi2);
#line 562 "slaln2.f"
	    bbnd = max(d__1,d__2);
#line 565 "slaln2.f"
	    if (bbnd > 1. && u22abs < 1.) {
#line 566 "slaln2.f"
		if (bbnd >= bignum * u22abs) {
#line 567 "slaln2.f"
		    *scale = 1. / bbnd;
#line 568 "slaln2.f"
		    br1 = *scale * br1;
#line 569 "slaln2.f"
		    bi1 = *scale * bi1;
#line 570 "slaln2.f"
		    br2 = *scale * br2;
#line 571 "slaln2.f"
		    bi2 = *scale * bi2;
#line 572 "slaln2.f"
		}
#line 573 "slaln2.f"
	    }

#line 575 "slaln2.f"
	    sladiv_(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
#line 576 "slaln2.f"
	    xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
#line 577 "slaln2.f"
	    xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
#line 578 "slaln2.f"
	    if (cswap[icmax - 1]) {
#line 579 "slaln2.f"
		x[x_dim1 + 1] = xr2;
#line 580 "slaln2.f"
		x[x_dim1 + 2] = xr1;
#line 581 "slaln2.f"
		x[(x_dim1 << 1) + 1] = xi2;
#line 582 "slaln2.f"
		x[(x_dim1 << 1) + 2] = xi1;
#line 583 "slaln2.f"
	    } else {
#line 584 "slaln2.f"
		x[x_dim1 + 1] = xr1;
#line 585 "slaln2.f"
		x[x_dim1 + 2] = xr2;
#line 586 "slaln2.f"
		x[(x_dim1 << 1) + 1] = xi1;
#line 587 "slaln2.f"
		x[(x_dim1 << 1) + 2] = xi2;
#line 588 "slaln2.f"
	    }
/* Computing MAX */
#line 589 "slaln2.f"
	    d__1 = abs(xr1) + abs(xi1), d__2 = abs(xr2) + abs(xi2);
#line 589 "slaln2.f"
	    *xnorm = max(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

#line 593 "slaln2.f"
	    if (*xnorm > 1. && cmax > 1.) {
#line 594 "slaln2.f"
		if (*xnorm > bignum / cmax) {
#line 595 "slaln2.f"
		    temp = cmax / bignum;
#line 596 "slaln2.f"
		    x[x_dim1 + 1] = temp * x[x_dim1 + 1];
#line 597 "slaln2.f"
		    x[x_dim1 + 2] = temp * x[x_dim1 + 2];
#line 598 "slaln2.f"
		    x[(x_dim1 << 1) + 1] = temp * x[(x_dim1 << 1) + 1];
#line 599 "slaln2.f"
		    x[(x_dim1 << 1) + 2] = temp * x[(x_dim1 << 1) + 2];
#line 600 "slaln2.f"
		    *xnorm = temp * *xnorm;
#line 601 "slaln2.f"
		    *scale = temp * *scale;
#line 602 "slaln2.f"
		}
#line 603 "slaln2.f"
	    }
#line 604 "slaln2.f"
	}
#line 605 "slaln2.f"
    }

#line 607 "slaln2.f"
    return 0;

/*     End of SLALN2 */

} /* slaln2_ */

#undef crv
#undef civ
#undef cr
#undef ci


