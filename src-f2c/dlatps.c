#line 1 "dlatps.f"
/* dlatps.f -- translated by f2c (version 20100827).
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

#line 1 "dlatps.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b DLATPS solves a triangular system of equations with the matrix held in packed storage. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLATPS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatps.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatps.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatps.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, */
/*                          CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), CNORM( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLATPS solves one of the triangular systems */
/* > */
/* >    A *x = s*b  or  A**T*x = s*b */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular matrix stored in packed form.  Here A**T denotes the */
/* > transpose of A, x and b are n-element vectors, and s is a scaling */
/* > factor, usually less than or equal to 1, chosen so that the */
/* > components of x will be less than the overflow threshold.  If the */
/* > unscaled problem will not cause overflow, the Level 2 BLAS routine */
/* > DTPSV is called. If the matrix A is singular (A(j,j) = 0 for some j), */
/* > then s is set to 0 and a non-trivial solution to A*x = 0 is returned. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the matrix A is upper or lower triangular. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the operation applied to A. */
/* >          = 'N':  Solve A * x = s*b  (No transpose) */
/* >          = 'T':  Solve A**T* x = s*b  (Transpose) */
/* >          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          Specifies whether or not the matrix A is unit triangular. */
/* >          = 'N':  Non-unit triangular */
/* >          = 'U':  Unit triangular */
/* > \endverbatim */
/* > */
/* > \param[in] NORMIN */
/* > \verbatim */
/* >          NORMIN is CHARACTER*1 */
/* >          Specifies whether CNORM has been set or not. */
/* >          = 'Y':  CNORM contains the column norms on entry */
/* >          = 'N':  CNORM is not set on entry.  On exit, the norms will */
/* >                  be computed and stored in CNORM. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the right hand side b of the triangular system. */
/* >          On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          The scaling factor s for the triangular system */
/* >             A * x = s*b  or  A**T* x = s*b. */
/* >          If SCALE = 0, the matrix A is singular or badly scaled, and */
/* >          the vector x is an exact or approximate solution to A*x = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CNORM */
/* > \verbatim */
/* >          CNORM is DOUBLE PRECISION array, dimension (N) */
/* > */
/* >          If NORMIN = 'Y', CNORM is an input argument and CNORM(j) */
/* >          contains the norm of the off-diagonal part of the j-th column */
/* >          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal */
/* >          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j) */
/* >          must be greater than or equal to the 1-norm. */
/* > */
/* >          If NORMIN = 'N', CNORM is an output argument and CNORM(j) */
/* >          returns the 1-norm of the offdiagonal part of the j-th column */
/* >          of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -k, the k-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, DTPSV */
/* >  is called, otherwise, specific code is used which checks for possible */
/* >  overflow or divide-by-zero at every operation. */
/* > */
/* >  A columnwise scheme is used for solving A*x = b.  The basic algorithm */
/* >  if A is lower triangular is */
/* > */
/* >       x[1:n] := b[1:n] */
/* >       for j = 1, ..., n */
/* >            x(j) := x(j) / A(j,j) */
/* >            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j] */
/* >       end */
/* > */
/* >  Define bounds on the components of x after j iterations of the loop: */
/* >     M(j) = bound on x[1:j] */
/* >     G(j) = bound on x[j+1:n] */
/* >  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}. */
/* > */
/* >  Then for iteration j+1 we have */
/* >     M(j+1) <= G(j) / | A(j+1,j+1) | */
/* >     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] | */
/* >            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | ) */
/* > */
/* >  where CNORM(j+1) is greater than or equal to the infinity-norm of */
/* >  column j+1 of A, not counting the diagonal.  Hence */
/* > */
/* >     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | ) */
/* >                  1<=i<=j */
/* >  and */
/* > */
/* >     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| ) */
/* >                                   1<=i< j */
/* > */
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTPSV if the */
/* >  reciprocal of the largest M(j), j=1,..,n, is larger than */
/* >  max(underflow, 1/overflow). */
/* > */
/* >  The bound on x(j) is also used to determine when a step in the */
/* >  columnwise method can be performed without fear of overflow.  If */
/* >  the computed bound is greater than a large constant, x is scaled to */
/* >  prevent overflow, but if the bound overflows, x is set to 0, x(j) to */
/* >  1, and scale to 0, and a non-trivial solution to A*x = 0 is found. */
/* > */
/* >  Similarly, a row-wise scheme is used to solve A**T*x = b.  The basic */
/* >  algorithm for A upper triangular is */
/* > */
/* >       for j = 1, ..., n */
/* >            x(j) := ( b(j) - A[1:j-1,j]**T * x[1:j-1] ) / A(j,j) */
/* >       end */
/* > */
/* >  We simultaneously compute two bounds */
/* >       G(j) = bound on ( b(i) - A[1:i-1,i]**T * x[1:i-1] ), 1<=i<=j */
/* >       M(j) = bound on x(i), 1<=i<=j */
/* > */
/* >  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we */
/* >  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1. */
/* >  Then the bound on x(j) is */
/* > */
/* >       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) | */
/* > */
/* >            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| ) */
/* >                      1<=i<=j */
/* > */
/* >  and we can safely call DTPSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *ap, doublereal *x, doublereal *scale, 
	doublereal *cnorm, integer *info, ftnlen uplo_len, ftnlen trans_len, 
	ftnlen diag_len, ftnlen normin_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, ip;
    static doublereal xj, rec, tjj;
    static integer jinc, jlen;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax, tjjs, xmax, grow, sumj;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal tscal, uscal;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer jlast;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int dtpsv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical notran;
    static integer jfirst;
    static doublereal smlnum;
    static logical nounit;


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
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 272 "dlatps.f"
    /* Parameter adjustments */
#line 272 "dlatps.f"
    --cnorm;
#line 272 "dlatps.f"
    --x;
#line 272 "dlatps.f"
    --ap;
#line 272 "dlatps.f"

#line 272 "dlatps.f"
    /* Function Body */
#line 272 "dlatps.f"
    *info = 0;
#line 273 "dlatps.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 274 "dlatps.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 275 "dlatps.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 279 "dlatps.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 280 "dlatps.f"
	*info = -1;
#line 281 "dlatps.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 283 "dlatps.f"
	*info = -2;
#line 284 "dlatps.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 285 "dlatps.f"
	*info = -3;
#line 286 "dlatps.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 288 "dlatps.f"
	*info = -4;
#line 289 "dlatps.f"
    } else if (*n < 0) {
#line 290 "dlatps.f"
	*info = -5;
#line 291 "dlatps.f"
    }
#line 292 "dlatps.f"
    if (*info != 0) {
#line 293 "dlatps.f"
	i__1 = -(*info);
#line 293 "dlatps.f"
	xerbla_("DLATPS", &i__1, (ftnlen)6);
#line 294 "dlatps.f"
	return 0;
#line 295 "dlatps.f"
    }

/*     Quick return if possible */

#line 299 "dlatps.f"
    if (*n == 0) {
#line 299 "dlatps.f"
	return 0;
#line 299 "dlatps.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 304 "dlatps.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 305 "dlatps.f"
    bignum = 1. / smlnum;
#line 306 "dlatps.f"
    *scale = 1.;

#line 308 "dlatps.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 312 "dlatps.f"
	if (upper) {

/*           A is upper triangular. */

#line 316 "dlatps.f"
	    ip = 1;
#line 317 "dlatps.f"
	    i__1 = *n;
#line 317 "dlatps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 318 "dlatps.f"
		i__2 = j - 1;
#line 318 "dlatps.f"
		cnorm[j] = dasum_(&i__2, &ap[ip], &c__1);
#line 319 "dlatps.f"
		ip += j;
#line 320 "dlatps.f"
/* L10: */
#line 320 "dlatps.f"
	    }
#line 321 "dlatps.f"
	} else {

/*           A is lower triangular. */

#line 325 "dlatps.f"
	    ip = 1;
#line 326 "dlatps.f"
	    i__1 = *n - 1;
#line 326 "dlatps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 327 "dlatps.f"
		i__2 = *n - j;
#line 327 "dlatps.f"
		cnorm[j] = dasum_(&i__2, &ap[ip + 1], &c__1);
#line 328 "dlatps.f"
		ip = ip + *n - j + 1;
#line 329 "dlatps.f"
/* L20: */
#line 329 "dlatps.f"
	    }
#line 330 "dlatps.f"
	    cnorm[*n] = 0.;
#line 331 "dlatps.f"
	}
#line 332 "dlatps.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM. */

#line 337 "dlatps.f"
    imax = idamax_(n, &cnorm[1], &c__1);
#line 338 "dlatps.f"
    tmax = cnorm[imax];
#line 339 "dlatps.f"
    if (tmax <= bignum) {
#line 340 "dlatps.f"
	tscal = 1.;
#line 341 "dlatps.f"
    } else {
#line 342 "dlatps.f"
	tscal = 1. / (smlnum * tmax);
#line 343 "dlatps.f"
	dscal_(n, &tscal, &cnorm[1], &c__1);
#line 344 "dlatps.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine DTPSV can be used. */

#line 349 "dlatps.f"
    j = idamax_(n, &x[1], &c__1);
#line 350 "dlatps.f"
    xmax = (d__1 = x[j], abs(d__1));
#line 351 "dlatps.f"
    xbnd = xmax;
#line 352 "dlatps.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 356 "dlatps.f"
	if (upper) {
#line 357 "dlatps.f"
	    jfirst = *n;
#line 358 "dlatps.f"
	    jlast = 1;
#line 359 "dlatps.f"
	    jinc = -1;
#line 360 "dlatps.f"
	} else {
#line 361 "dlatps.f"
	    jfirst = 1;
#line 362 "dlatps.f"
	    jlast = *n;
#line 363 "dlatps.f"
	    jinc = 1;
#line 364 "dlatps.f"
	}

#line 366 "dlatps.f"
	if (tscal != 1.) {
#line 367 "dlatps.f"
	    grow = 0.;
#line 368 "dlatps.f"
	    goto L50;
#line 369 "dlatps.f"
	}

#line 371 "dlatps.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 378 "dlatps.f"
	    grow = 1. / max(xbnd,smlnum);
#line 379 "dlatps.f"
	    xbnd = grow;
#line 380 "dlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 381 "dlatps.f"
	    jlen = *n;
#line 382 "dlatps.f"
	    i__1 = jlast;
#line 382 "dlatps.f"
	    i__2 = jinc;
#line 382 "dlatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 386 "dlatps.f"
		if (grow <= smlnum) {
#line 386 "dlatps.f"
		    goto L50;
#line 386 "dlatps.f"
		}

/*              M(j) = G(j-1) / abs(A(j,j)) */

#line 391 "dlatps.f"
		tjj = (d__1 = ap[ip], abs(d__1));
/* Computing MIN */
#line 392 "dlatps.f"
		d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 392 "dlatps.f"
		xbnd = min(d__1,d__2);
#line 393 "dlatps.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 397 "dlatps.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 398 "dlatps.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 402 "dlatps.f"
		    grow = 0.;
#line 403 "dlatps.f"
		}
#line 404 "dlatps.f"
		ip += jinc * jlen;
#line 405 "dlatps.f"
		--jlen;
#line 406 "dlatps.f"
/* L30: */
#line 406 "dlatps.f"
	    }
#line 407 "dlatps.f"
	    grow = xbnd;
#line 408 "dlatps.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 414 "dlatps.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 414 "dlatps.f"
	    grow = min(d__1,d__2);
#line 415 "dlatps.f"
	    i__2 = jlast;
#line 415 "dlatps.f"
	    i__1 = jinc;
#line 415 "dlatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 419 "dlatps.f"
		if (grow <= smlnum) {
#line 419 "dlatps.f"
		    goto L50;
#line 419 "dlatps.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 424 "dlatps.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 425 "dlatps.f"
/* L40: */
#line 425 "dlatps.f"
	    }
#line 426 "dlatps.f"
	}
#line 427 "dlatps.f"
L50:

#line 429 "dlatps.f"
	;
#line 429 "dlatps.f"
    } else {

/*        Compute the growth in A**T * x = b. */

#line 433 "dlatps.f"
	if (upper) {
#line 434 "dlatps.f"
	    jfirst = 1;
#line 435 "dlatps.f"
	    jlast = *n;
#line 436 "dlatps.f"
	    jinc = 1;
#line 437 "dlatps.f"
	} else {
#line 438 "dlatps.f"
	    jfirst = *n;
#line 439 "dlatps.f"
	    jlast = 1;
#line 440 "dlatps.f"
	    jinc = -1;
#line 441 "dlatps.f"
	}

#line 443 "dlatps.f"
	if (tscal != 1.) {
#line 444 "dlatps.f"
	    grow = 0.;
#line 445 "dlatps.f"
	    goto L80;
#line 446 "dlatps.f"
	}

#line 448 "dlatps.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 455 "dlatps.f"
	    grow = 1. / max(xbnd,smlnum);
#line 456 "dlatps.f"
	    xbnd = grow;
#line 457 "dlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 458 "dlatps.f"
	    jlen = 1;
#line 459 "dlatps.f"
	    i__1 = jlast;
#line 459 "dlatps.f"
	    i__2 = jinc;
#line 459 "dlatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 463 "dlatps.f"
		if (grow <= smlnum) {
#line 463 "dlatps.f"
		    goto L80;
#line 463 "dlatps.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 468 "dlatps.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 469 "dlatps.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 469 "dlatps.f"
		grow = min(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 473 "dlatps.f"
		tjj = (d__1 = ap[ip], abs(d__1));
#line 474 "dlatps.f"
		if (xj > tjj) {
#line 474 "dlatps.f"
		    xbnd *= tjj / xj;
#line 474 "dlatps.f"
		}
#line 476 "dlatps.f"
		++jlen;
#line 477 "dlatps.f"
		ip += jinc * jlen;
#line 478 "dlatps.f"
/* L60: */
#line 478 "dlatps.f"
	    }
#line 479 "dlatps.f"
	    grow = min(grow,xbnd);
#line 480 "dlatps.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 486 "dlatps.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 486 "dlatps.f"
	    grow = min(d__1,d__2);
#line 487 "dlatps.f"
	    i__2 = jlast;
#line 487 "dlatps.f"
	    i__1 = jinc;
#line 487 "dlatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 491 "dlatps.f"
		if (grow <= smlnum) {
#line 491 "dlatps.f"
		    goto L80;
#line 491 "dlatps.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 496 "dlatps.f"
		xj = cnorm[j] + 1.;
#line 497 "dlatps.f"
		grow /= xj;
#line 498 "dlatps.f"
/* L70: */
#line 498 "dlatps.f"
	    }
#line 499 "dlatps.f"
	}
#line 500 "dlatps.f"
L80:
#line 501 "dlatps.f"
	;
#line 501 "dlatps.f"
    }

#line 503 "dlatps.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 508 "dlatps.f"
	dtpsv_(uplo, trans, diag, n, &ap[1], &x[1], &c__1, (ftnlen)1, (ftnlen)
		1, (ftnlen)1);
#line 509 "dlatps.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 513 "dlatps.f"
	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 518 "dlatps.f"
	    *scale = bignum / xmax;
#line 519 "dlatps.f"
	    dscal_(n, scale, &x[1], &c__1);
#line 520 "dlatps.f"
	    xmax = bignum;
#line 521 "dlatps.f"
	}

#line 523 "dlatps.f"
	if (notran) {

/*           Solve A * x = b */

#line 527 "dlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 528 "dlatps.f"
	    i__1 = jlast;
#line 528 "dlatps.f"
	    i__2 = jinc;
#line 528 "dlatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 532 "dlatps.f"
		xj = (d__1 = x[j], abs(d__1));
#line 533 "dlatps.f"
		if (nounit) {
#line 534 "dlatps.f"
		    tjjs = ap[ip] * tscal;
#line 535 "dlatps.f"
		} else {
#line 536 "dlatps.f"
		    tjjs = tscal;
#line 537 "dlatps.f"
		    if (tscal == 1.) {
#line 537 "dlatps.f"
			goto L100;
#line 537 "dlatps.f"
		    }
#line 539 "dlatps.f"
		}
#line 540 "dlatps.f"
		tjj = abs(tjjs);
#line 541 "dlatps.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 545 "dlatps.f"
		    if (tjj < 1.) {
#line 546 "dlatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 550 "dlatps.f"
			    rec = 1. / xj;
#line 551 "dlatps.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 552 "dlatps.f"
			    *scale *= rec;
#line 553 "dlatps.f"
			    xmax *= rec;
#line 554 "dlatps.f"
			}
#line 555 "dlatps.f"
		    }
#line 556 "dlatps.f"
		    x[j] /= tjjs;
#line 557 "dlatps.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 558 "dlatps.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 562 "dlatps.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 567 "dlatps.f"
			rec = tjj * bignum / xj;
#line 568 "dlatps.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 573 "dlatps.f"
			    rec /= cnorm[j];
#line 574 "dlatps.f"
			}
#line 575 "dlatps.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 576 "dlatps.f"
			*scale *= rec;
#line 577 "dlatps.f"
			xmax *= rec;
#line 578 "dlatps.f"
		    }
#line 579 "dlatps.f"
		    x[j] /= tjjs;
#line 580 "dlatps.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 581 "dlatps.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 586 "dlatps.f"
		    i__3 = *n;
#line 586 "dlatps.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 587 "dlatps.f"
			x[i__] = 0.;
#line 588 "dlatps.f"
/* L90: */
#line 588 "dlatps.f"
		    }
#line 589 "dlatps.f"
		    x[j] = 1.;
#line 590 "dlatps.f"
		    xj = 1.;
#line 591 "dlatps.f"
		    *scale = 0.;
#line 592 "dlatps.f"
		    xmax = 0.;
#line 593 "dlatps.f"
		}
#line 594 "dlatps.f"
L100:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 599 "dlatps.f"
		if (xj > 1.) {
#line 600 "dlatps.f"
		    rec = 1. / xj;
#line 601 "dlatps.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 605 "dlatps.f"
			rec *= .5;
#line 606 "dlatps.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 607 "dlatps.f"
			*scale *= rec;
#line 608 "dlatps.f"
		    }
#line 609 "dlatps.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 613 "dlatps.f"
		    dscal_(n, &c_b36, &x[1], &c__1);
#line 614 "dlatps.f"
		    *scale *= .5;
#line 615 "dlatps.f"
		}

#line 617 "dlatps.f"
		if (upper) {
#line 618 "dlatps.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

#line 623 "dlatps.f"
			i__3 = j - 1;
#line 623 "dlatps.f"
			d__1 = -x[j] * tscal;
#line 623 "dlatps.f"
			daxpy_(&i__3, &d__1, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 625 "dlatps.f"
			i__3 = j - 1;
#line 625 "dlatps.f"
			i__ = idamax_(&i__3, &x[1], &c__1);
#line 626 "dlatps.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 627 "dlatps.f"
		    }
#line 628 "dlatps.f"
		    ip -= j;
#line 629 "dlatps.f"
		} else {
#line 630 "dlatps.f"
		    if (j < *n) {

/*                    Compute the update */
/*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

#line 635 "dlatps.f"
			i__3 = *n - j;
#line 635 "dlatps.f"
			d__1 = -x[j] * tscal;
#line 635 "dlatps.f"
			daxpy_(&i__3, &d__1, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 637 "dlatps.f"
			i__3 = *n - j;
#line 637 "dlatps.f"
			i__ = j + idamax_(&i__3, &x[j + 1], &c__1);
#line 638 "dlatps.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 639 "dlatps.f"
		    }
#line 640 "dlatps.f"
		    ip = ip + *n - j + 1;
#line 641 "dlatps.f"
		}
#line 642 "dlatps.f"
/* L110: */
#line 642 "dlatps.f"
	    }

#line 644 "dlatps.f"
	} else {

/*           Solve A**T * x = b */

#line 648 "dlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 649 "dlatps.f"
	    jlen = 1;
#line 650 "dlatps.f"
	    i__2 = jlast;
#line 650 "dlatps.f"
	    i__1 = jinc;
#line 650 "dlatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 655 "dlatps.f"
		xj = (d__1 = x[j], abs(d__1));
#line 656 "dlatps.f"
		uscal = tscal;
#line 657 "dlatps.f"
		rec = 1. / max(xmax,1.);
#line 658 "dlatps.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 662 "dlatps.f"
		    rec *= .5;
#line 663 "dlatps.f"
		    if (nounit) {
#line 664 "dlatps.f"
			tjjs = ap[ip] * tscal;
#line 665 "dlatps.f"
		    } else {
#line 666 "dlatps.f"
			tjjs = tscal;
#line 667 "dlatps.f"
		    }
#line 668 "dlatps.f"
		    tjj = abs(tjjs);
#line 669 "dlatps.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 673 "dlatps.f"
			d__1 = 1., d__2 = rec * tjj;
#line 673 "dlatps.f"
			rec = min(d__1,d__2);
#line 674 "dlatps.f"
			uscal /= tjjs;
#line 675 "dlatps.f"
		    }
#line 676 "dlatps.f"
		    if (rec < 1.) {
#line 677 "dlatps.f"
			dscal_(n, &rec, &x[1], &c__1);
#line 678 "dlatps.f"
			*scale *= rec;
#line 679 "dlatps.f"
			xmax *= rec;
#line 680 "dlatps.f"
		    }
#line 681 "dlatps.f"
		}

#line 683 "dlatps.f"
		sumj = 0.;
#line 684 "dlatps.f"
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call DDOT to perform the dot product. */

#line 689 "dlatps.f"
		    if (upper) {
#line 690 "dlatps.f"
			i__3 = j - 1;
#line 690 "dlatps.f"
			sumj = ddot_(&i__3, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 691 "dlatps.f"
		    } else if (j < *n) {
#line 692 "dlatps.f"
			i__3 = *n - j;
#line 692 "dlatps.f"
			sumj = ddot_(&i__3, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 693 "dlatps.f"
		    }
#line 694 "dlatps.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 698 "dlatps.f"
		    if (upper) {
#line 699 "dlatps.f"
			i__3 = j - 1;
#line 699 "dlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 700 "dlatps.f"
			    sumj += ap[ip - j + i__] * uscal * x[i__];
#line 701 "dlatps.f"
/* L120: */
#line 701 "dlatps.f"
			}
#line 702 "dlatps.f"
		    } else if (j < *n) {
#line 703 "dlatps.f"
			i__3 = *n - j;
#line 703 "dlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 704 "dlatps.f"
			    sumj += ap[ip + i__] * uscal * x[j + i__];
#line 705 "dlatps.f"
/* L130: */
#line 705 "dlatps.f"
			}
#line 706 "dlatps.f"
		    }
#line 707 "dlatps.f"
		}

#line 709 "dlatps.f"
		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 714 "dlatps.f"
		    x[j] -= sumj;
#line 715 "dlatps.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 716 "dlatps.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 720 "dlatps.f"
			tjjs = ap[ip] * tscal;
#line 721 "dlatps.f"
		    } else {
#line 722 "dlatps.f"
			tjjs = tscal;
#line 723 "dlatps.f"
			if (tscal == 1.) {
#line 723 "dlatps.f"
			    goto L150;
#line 723 "dlatps.f"
			}
#line 725 "dlatps.f"
		    }
#line 726 "dlatps.f"
		    tjj = abs(tjjs);
#line 727 "dlatps.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 731 "dlatps.f"
			if (tjj < 1.) {
#line 732 "dlatps.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 736 "dlatps.f"
				rec = 1. / xj;
#line 737 "dlatps.f"
				dscal_(n, &rec, &x[1], &c__1);
#line 738 "dlatps.f"
				*scale *= rec;
#line 739 "dlatps.f"
				xmax *= rec;
#line 740 "dlatps.f"
			    }
#line 741 "dlatps.f"
			}
#line 742 "dlatps.f"
			x[j] /= tjjs;
#line 743 "dlatps.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 747 "dlatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 751 "dlatps.f"
			    rec = tjj * bignum / xj;
#line 752 "dlatps.f"
			    dscal_(n, &rec, &x[1], &c__1);
#line 753 "dlatps.f"
			    *scale *= rec;
#line 754 "dlatps.f"
			    xmax *= rec;
#line 755 "dlatps.f"
			}
#line 756 "dlatps.f"
			x[j] /= tjjs;
#line 757 "dlatps.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0, and compute a solution to A**T*x = 0. */

#line 762 "dlatps.f"
			i__3 = *n;
#line 762 "dlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 763 "dlatps.f"
			    x[i__] = 0.;
#line 764 "dlatps.f"
/* L140: */
#line 764 "dlatps.f"
			}
#line 765 "dlatps.f"
			x[j] = 1.;
#line 766 "dlatps.f"
			*scale = 0.;
#line 767 "dlatps.f"
			xmax = 0.;
#line 768 "dlatps.f"
		    }
#line 769 "dlatps.f"
L150:
#line 770 "dlatps.f"
		    ;
#line 770 "dlatps.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 775 "dlatps.f"
		    x[j] = x[j] / tjjs - sumj;
#line 776 "dlatps.f"
		}
/* Computing MAX */
#line 777 "dlatps.f"
		d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
#line 777 "dlatps.f"
		xmax = max(d__2,d__3);
#line 778 "dlatps.f"
		++jlen;
#line 779 "dlatps.f"
		ip += jinc * jlen;
#line 780 "dlatps.f"
/* L160: */
#line 780 "dlatps.f"
	    }
#line 781 "dlatps.f"
	}
#line 782 "dlatps.f"
	*scale /= tscal;
#line 783 "dlatps.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 787 "dlatps.f"
    if (tscal != 1.) {
#line 788 "dlatps.f"
	d__1 = 1. / tscal;
#line 788 "dlatps.f"
	dscal_(n, &d__1, &cnorm[1], &c__1);
#line 789 "dlatps.f"
    }

#line 791 "dlatps.f"
    return 0;

/*     End of DLATPS */

} /* dlatps_ */

