#line 1 "slasd8.f"
/* slasd8.f -- translated by f2c (version 20100827).
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

#line 1 "slasd8.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b8 = 1.;

/* > \brief \b SLASD8 finds the square roots of the roots of the secular equation, and stores, for each elemen
t in D, the distance to its two nearest poles. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASD8 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd8.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd8.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd8.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, */
/*                          DSIGMA, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ICOMPQ, INFO, K, LDDIFR */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), DIFL( * ), DIFR( LDDIFR, * ), */
/*      $                   DSIGMA( * ), VF( * ), VL( * ), WORK( * ), */
/*      $                   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASD8 finds the square roots of the roots of the secular equation, */
/* > as defined by the values in DSIGMA and Z. It makes the appropriate */
/* > calls to SLASD4, and stores, for each  element in D, the distance */
/* > to its two nearest poles (elements in DSIGMA). It also updates */
/* > the arrays VF and VL, the first and last components of all the */
/* > right singular vectors of the original bidiagonal matrix. */
/* > */
/* > SLASD8 is called from SLASD6. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ICOMPQ */
/* > \verbatim */
/* >          ICOMPQ is INTEGER */
/* >          Specifies whether singular vectors are to be computed in */
/* >          factored form in the calling routine: */
/* >          = 0: Compute singular values only. */
/* >          = 1: Compute singular vectors in factored form as well. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of terms in the rational function to be solved */
/* >          by SLASD4.  K >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is REAL array, dimension ( K ) */
/* >          On output, D contains the updated singular values. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension ( K ) */
/* >          On entry, the first K elements of this array contain the */
/* >          components of the deflation-adjusted updating row vector. */
/* >          On exit, Z is updated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VF */
/* > \verbatim */
/* >          VF is REAL array, dimension ( K ) */
/* >          On entry, VF contains  information passed through DBEDE8. */
/* >          On exit, VF contains the first K components of the first */
/* >          components of all right singular vectors of the bidiagonal */
/* >          matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is REAL array, dimension ( K ) */
/* >          On entry, VL contains  information passed through DBEDE8. */
/* >          On exit, VL contains the first K components of the last */
/* >          components of all right singular vectors of the bidiagonal */
/* >          matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] DIFL */
/* > \verbatim */
/* >          DIFL is REAL array, dimension ( K ) */
/* >          On exit, DIFL(I) = D(I) - DSIGMA(I). */
/* > \endverbatim */
/* > */
/* > \param[out] DIFR */
/* > \verbatim */
/* >          DIFR is REAL array, */
/* >                   dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and */
/* >                   dimension ( K ) if ICOMPQ = 0. */
/* >          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not */
/* >          defined and will not be referenced. */
/* > */
/* >          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the */
/* >          normalizing factors for the right singular vector matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDDIFR */
/* > \verbatim */
/* >          LDDIFR is INTEGER */
/* >          The leading dimension of DIFR, must be at least K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DSIGMA */
/* > \verbatim */
/* >          DSIGMA is REAL array, dimension ( K ) */
/* >          On entry, the first K elements of this array contain the old */
/* >          roots of the deflated updating problem.  These are the poles */
/* >          of the secular equation. */
/* >          On exit, the elements of DSIGMA may be very slightly altered */
/* >          in value. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension at least 3 * K */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = 1, a singular value did not converge */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slasd8_(integer *icompq, integer *k, doublereal *d__, 
	doublereal *z__, doublereal *vf, doublereal *vl, doublereal *difl, 
	doublereal *difr, integer *lddifr, doublereal *dsigma, doublereal *
	work, integer *info)
{
    /* System generated locals */
    integer difr_dim1, difr_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j;
    static doublereal dj, rho;
    static integer iwk1, iwk2, iwk3;
    static doublereal temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer iwk2i, iwk3i;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal diflj, difrj, dsigj;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal slamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int slasd4_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);
    static doublereal dsigjp;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), slaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen);


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
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 207 "slasd8.f"
    /* Parameter adjustments */
#line 207 "slasd8.f"
    --d__;
#line 207 "slasd8.f"
    --z__;
#line 207 "slasd8.f"
    --vf;
#line 207 "slasd8.f"
    --vl;
#line 207 "slasd8.f"
    --difl;
#line 207 "slasd8.f"
    difr_dim1 = *lddifr;
#line 207 "slasd8.f"
    difr_offset = 1 + difr_dim1;
#line 207 "slasd8.f"
    difr -= difr_offset;
#line 207 "slasd8.f"
    --dsigma;
#line 207 "slasd8.f"
    --work;
#line 207 "slasd8.f"

#line 207 "slasd8.f"
    /* Function Body */
#line 207 "slasd8.f"
    *info = 0;

#line 209 "slasd8.f"
    if (*icompq < 0 || *icompq > 1) {
#line 210 "slasd8.f"
	*info = -1;
#line 211 "slasd8.f"
    } else if (*k < 1) {
#line 212 "slasd8.f"
	*info = -2;
#line 213 "slasd8.f"
    } else if (*lddifr < *k) {
#line 214 "slasd8.f"
	*info = -9;
#line 215 "slasd8.f"
    }
#line 216 "slasd8.f"
    if (*info != 0) {
#line 217 "slasd8.f"
	i__1 = -(*info);
#line 217 "slasd8.f"
	xerbla_("SLASD8", &i__1, (ftnlen)6);
#line 218 "slasd8.f"
	return 0;
#line 219 "slasd8.f"
    }

/*     Quick return if possible */

#line 223 "slasd8.f"
    if (*k == 1) {
#line 224 "slasd8.f"
	d__[1] = abs(z__[1]);
#line 225 "slasd8.f"
	difl[1] = d__[1];
#line 226 "slasd8.f"
	if (*icompq == 1) {
#line 227 "slasd8.f"
	    difl[2] = 1.;
#line 228 "slasd8.f"
	    difr[(difr_dim1 << 1) + 1] = 1.;
#line 229 "slasd8.f"
	}
#line 230 "slasd8.f"
	return 0;
#line 231 "slasd8.f"
    }

/*     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can */
/*     be computed with high relative accuracy (barring over/underflow). */
/*     This is a problem on machines without a guard digit in */
/*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2). */
/*     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I), */
/*     which on any of these machines zeros out the bottommost */
/*     bit of DSIGMA(I) if it is 1; this makes the subsequent */
/*     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation */
/*     occurs. On binary machines with a guard digit (almost all */
/*     machines) it does not change DSIGMA(I) at all. On hexadecimal */
/*     and decimal machines with a guard digit, it slightly */
/*     changes the bottommost bits of DSIGMA(I). It does not account */
/*     for hexadecimal or decimal machines without guard digits */
/*     (we know of none). We use a subroutine call to compute */
/*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating */
/*     this code. */

#line 250 "slasd8.f"
    i__1 = *k;
#line 250 "slasd8.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "slasd8.f"
	dsigma[i__] = slamc3_(&dsigma[i__], &dsigma[i__]) - dsigma[i__];
#line 252 "slasd8.f"
/* L10: */
#line 252 "slasd8.f"
    }

/*     Book keeping. */

#line 256 "slasd8.f"
    iwk1 = 1;
#line 257 "slasd8.f"
    iwk2 = iwk1 + *k;
#line 258 "slasd8.f"
    iwk3 = iwk2 + *k;
#line 259 "slasd8.f"
    iwk2i = iwk2 - 1;
#line 260 "slasd8.f"
    iwk3i = iwk3 - 1;

/*     Normalize Z. */

#line 264 "slasd8.f"
    rho = snrm2_(k, &z__[1], &c__1);
#line 265 "slasd8.f"
    slascl_("G", &c__0, &c__0, &rho, &c_b8, k, &c__1, &z__[1], k, info, (
	    ftnlen)1);
#line 266 "slasd8.f"
    rho *= rho;

/*     Initialize WORK(IWK3). */

#line 270 "slasd8.f"
    slaset_("A", k, &c__1, &c_b8, &c_b8, &work[iwk3], k, (ftnlen)1);

/*     Compute the updated singular values, the arrays DIFL, DIFR, */
/*     and the updated Z. */

#line 275 "slasd8.f"
    i__1 = *k;
#line 275 "slasd8.f"
    for (j = 1; j <= i__1; ++j) {
#line 276 "slasd8.f"
	slasd4_(k, &j, &dsigma[1], &z__[1], &work[iwk1], &rho, &d__[j], &work[
		iwk2], info);

/*        If the root finder fails, report the convergence failure. */

#line 281 "slasd8.f"
	if (*info != 0) {
#line 282 "slasd8.f"
	    return 0;
#line 283 "slasd8.f"
	}
#line 284 "slasd8.f"
	work[iwk3i + j] = work[iwk3i + j] * work[j] * work[iwk2i + j];
#line 285 "slasd8.f"
	difl[j] = -work[j];
#line 286 "slasd8.f"
	difr[j + difr_dim1] = -work[j + 1];
#line 287 "slasd8.f"
	i__2 = j - 1;
#line 287 "slasd8.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 288 "slasd8.f"
	    work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + 
		    i__] / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[
		    j]);
#line 292 "slasd8.f"
/* L20: */
#line 292 "slasd8.f"
	}
#line 293 "slasd8.f"
	i__2 = *k;
#line 293 "slasd8.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 294 "slasd8.f"
	    work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + 
		    i__] / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[
		    j]);
#line 298 "slasd8.f"
/* L30: */
#line 298 "slasd8.f"
	}
#line 299 "slasd8.f"
/* L40: */
#line 299 "slasd8.f"
    }

/*     Compute updated Z. */

#line 303 "slasd8.f"
    i__1 = *k;
#line 303 "slasd8.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 304 "slasd8.f"
	d__2 = sqrt((d__1 = work[iwk3i + i__], abs(d__1)));
#line 304 "slasd8.f"
	z__[i__] = d_sign(&d__2, &z__[i__]);
#line 305 "slasd8.f"
/* L50: */
#line 305 "slasd8.f"
    }

/*     Update VF and VL. */

#line 309 "slasd8.f"
    i__1 = *k;
#line 309 "slasd8.f"
    for (j = 1; j <= i__1; ++j) {
#line 310 "slasd8.f"
	diflj = difl[j];
#line 311 "slasd8.f"
	dj = d__[j];
#line 312 "slasd8.f"
	dsigj = -dsigma[j];
#line 313 "slasd8.f"
	if (j < *k) {
#line 314 "slasd8.f"
	    difrj = -difr[j + difr_dim1];
#line 315 "slasd8.f"
	    dsigjp = -dsigma[j + 1];
#line 316 "slasd8.f"
	}
#line 317 "slasd8.f"
	work[j] = -z__[j] / diflj / (dsigma[j] + dj);
#line 318 "slasd8.f"
	i__2 = j - 1;
#line 318 "slasd8.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 319 "slasd8.f"
	    work[i__] = z__[i__] / (slamc3_(&dsigma[i__], &dsigj) - diflj) / (
		    dsigma[i__] + dj);
#line 321 "slasd8.f"
/* L60: */
#line 321 "slasd8.f"
	}
#line 322 "slasd8.f"
	i__2 = *k;
#line 322 "slasd8.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 323 "slasd8.f"
	    work[i__] = z__[i__] / (slamc3_(&dsigma[i__], &dsigjp) + difrj) / 
		    (dsigma[i__] + dj);
#line 325 "slasd8.f"
/* L70: */
#line 325 "slasd8.f"
	}
#line 326 "slasd8.f"
	temp = snrm2_(k, &work[1], &c__1);
#line 327 "slasd8.f"
	work[iwk2i + j] = sdot_(k, &work[1], &c__1, &vf[1], &c__1) / temp;
#line 328 "slasd8.f"
	work[iwk3i + j] = sdot_(k, &work[1], &c__1, &vl[1], &c__1) / temp;
#line 329 "slasd8.f"
	if (*icompq == 1) {
#line 330 "slasd8.f"
	    difr[j + (difr_dim1 << 1)] = temp;
#line 331 "slasd8.f"
	}
#line 332 "slasd8.f"
/* L80: */
#line 332 "slasd8.f"
    }

#line 334 "slasd8.f"
    scopy_(k, &work[iwk2], &c__1, &vf[1], &c__1);
#line 335 "slasd8.f"
    scopy_(k, &work[iwk3], &c__1, &vl[1], &c__1);

#line 337 "slasd8.f"
    return 0;

/*     End of SLASD8 */

} /* slasd8_ */

