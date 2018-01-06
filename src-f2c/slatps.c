#line 1 "slatps.f"
/* slatps.f -- translated by f2c (version 20100827).
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

#line 1 "slatps.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b SLATPS solves a triangular system of equations with the matrix held in packed storage. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLATPS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatps.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatps.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatps.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, */
/*                          CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ), CNORM( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLATPS solves one of the triangular systems */
/* > */
/* >    A *x = s*b  or  A**T*x = s*b */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular matrix stored in packed form.  Here A**T denotes the */
/* > transpose of A, x and b are n-element vectors, and s is a scaling */
/* > factor, usually less than or equal to 1, chosen so that the */
/* > components of x will be less than the overflow threshold.  If the */
/* > unscaled problem will not cause overflow, the Level 2 BLAS routine */
/* > STPSV is called. If the matrix A is singular (A(j,j) = 0 for some j), */
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
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (N) */
/* >          On entry, the right hand side b of the triangular system. */
/* >          On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >          The scaling factor s for the triangular system */
/* >             A * x = s*b  or  A**T* x = s*b. */
/* >          If SCALE = 0, the matrix A is singular or badly scaled, and */
/* >          the vector x is an exact or approximate solution to A*x = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CNORM */
/* > \verbatim */
/* >          CNORM is REAL array, dimension (N) */
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

/* > \ingroup realOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, STPSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine STPSV if the */
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
/* >  and we can safely call STPSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slatps_(char *uplo, char *trans, char *diag, char *
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
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax, tjjs;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal xmax, grow, sumj;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal tscal, uscal;
    static integer jlast;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), stpsv_(char *, char *, char *
	    , integer *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
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

#line 272 "slatps.f"
    /* Parameter adjustments */
#line 272 "slatps.f"
    --cnorm;
#line 272 "slatps.f"
    --x;
#line 272 "slatps.f"
    --ap;
#line 272 "slatps.f"

#line 272 "slatps.f"
    /* Function Body */
#line 272 "slatps.f"
    *info = 0;
#line 273 "slatps.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 274 "slatps.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 275 "slatps.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 279 "slatps.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 280 "slatps.f"
	*info = -1;
#line 281 "slatps.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 283 "slatps.f"
	*info = -2;
#line 284 "slatps.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 285 "slatps.f"
	*info = -3;
#line 286 "slatps.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 288 "slatps.f"
	*info = -4;
#line 289 "slatps.f"
    } else if (*n < 0) {
#line 290 "slatps.f"
	*info = -5;
#line 291 "slatps.f"
    }
#line 292 "slatps.f"
    if (*info != 0) {
#line 293 "slatps.f"
	i__1 = -(*info);
#line 293 "slatps.f"
	xerbla_("SLATPS", &i__1, (ftnlen)6);
#line 294 "slatps.f"
	return 0;
#line 295 "slatps.f"
    }

/*     Quick return if possible */

#line 299 "slatps.f"
    if (*n == 0) {
#line 299 "slatps.f"
	return 0;
#line 299 "slatps.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 304 "slatps.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 305 "slatps.f"
    bignum = 1. / smlnum;
#line 306 "slatps.f"
    *scale = 1.;

#line 308 "slatps.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 312 "slatps.f"
	if (upper) {

/*           A is upper triangular. */

#line 316 "slatps.f"
	    ip = 1;
#line 317 "slatps.f"
	    i__1 = *n;
#line 317 "slatps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 318 "slatps.f"
		i__2 = j - 1;
#line 318 "slatps.f"
		cnorm[j] = sasum_(&i__2, &ap[ip], &c__1);
#line 319 "slatps.f"
		ip += j;
#line 320 "slatps.f"
/* L10: */
#line 320 "slatps.f"
	    }
#line 321 "slatps.f"
	} else {

/*           A is lower triangular. */

#line 325 "slatps.f"
	    ip = 1;
#line 326 "slatps.f"
	    i__1 = *n - 1;
#line 326 "slatps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 327 "slatps.f"
		i__2 = *n - j;
#line 327 "slatps.f"
		cnorm[j] = sasum_(&i__2, &ap[ip + 1], &c__1);
#line 328 "slatps.f"
		ip = ip + *n - j + 1;
#line 329 "slatps.f"
/* L20: */
#line 329 "slatps.f"
	    }
#line 330 "slatps.f"
	    cnorm[*n] = 0.;
#line 331 "slatps.f"
	}
#line 332 "slatps.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM. */

#line 337 "slatps.f"
    imax = isamax_(n, &cnorm[1], &c__1);
#line 338 "slatps.f"
    tmax = cnorm[imax];
#line 339 "slatps.f"
    if (tmax <= bignum) {
#line 340 "slatps.f"
	tscal = 1.;
#line 341 "slatps.f"
    } else {
#line 342 "slatps.f"
	tscal = 1. / (smlnum * tmax);
#line 343 "slatps.f"
	sscal_(n, &tscal, &cnorm[1], &c__1);
#line 344 "slatps.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine STPSV can be used. */

#line 349 "slatps.f"
    j = isamax_(n, &x[1], &c__1);
#line 350 "slatps.f"
    xmax = (d__1 = x[j], abs(d__1));
#line 351 "slatps.f"
    xbnd = xmax;
#line 352 "slatps.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 356 "slatps.f"
	if (upper) {
#line 357 "slatps.f"
	    jfirst = *n;
#line 358 "slatps.f"
	    jlast = 1;
#line 359 "slatps.f"
	    jinc = -1;
#line 360 "slatps.f"
	} else {
#line 361 "slatps.f"
	    jfirst = 1;
#line 362 "slatps.f"
	    jlast = *n;
#line 363 "slatps.f"
	    jinc = 1;
#line 364 "slatps.f"
	}

#line 366 "slatps.f"
	if (tscal != 1.) {
#line 367 "slatps.f"
	    grow = 0.;
#line 368 "slatps.f"
	    goto L50;
#line 369 "slatps.f"
	}

#line 371 "slatps.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 378 "slatps.f"
	    grow = 1. / max(xbnd,smlnum);
#line 379 "slatps.f"
	    xbnd = grow;
#line 380 "slatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 381 "slatps.f"
	    jlen = *n;
#line 382 "slatps.f"
	    i__1 = jlast;
#line 382 "slatps.f"
	    i__2 = jinc;
#line 382 "slatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 386 "slatps.f"
		if (grow <= smlnum) {
#line 386 "slatps.f"
		    goto L50;
#line 386 "slatps.f"
		}

/*              M(j) = G(j-1) / abs(A(j,j)) */

#line 391 "slatps.f"
		tjj = (d__1 = ap[ip], abs(d__1));
/* Computing MIN */
#line 392 "slatps.f"
		d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 392 "slatps.f"
		xbnd = min(d__1,d__2);
#line 393 "slatps.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 397 "slatps.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 398 "slatps.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 402 "slatps.f"
		    grow = 0.;
#line 403 "slatps.f"
		}
#line 404 "slatps.f"
		ip += jinc * jlen;
#line 405 "slatps.f"
		--jlen;
#line 406 "slatps.f"
/* L30: */
#line 406 "slatps.f"
	    }
#line 407 "slatps.f"
	    grow = xbnd;
#line 408 "slatps.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 414 "slatps.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 414 "slatps.f"
	    grow = min(d__1,d__2);
#line 415 "slatps.f"
	    i__2 = jlast;
#line 415 "slatps.f"
	    i__1 = jinc;
#line 415 "slatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 419 "slatps.f"
		if (grow <= smlnum) {
#line 419 "slatps.f"
		    goto L50;
#line 419 "slatps.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 424 "slatps.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 425 "slatps.f"
/* L40: */
#line 425 "slatps.f"
	    }
#line 426 "slatps.f"
	}
#line 427 "slatps.f"
L50:

#line 429 "slatps.f"
	;
#line 429 "slatps.f"
    } else {

/*        Compute the growth in A**T * x = b. */

#line 433 "slatps.f"
	if (upper) {
#line 434 "slatps.f"
	    jfirst = 1;
#line 435 "slatps.f"
	    jlast = *n;
#line 436 "slatps.f"
	    jinc = 1;
#line 437 "slatps.f"
	} else {
#line 438 "slatps.f"
	    jfirst = *n;
#line 439 "slatps.f"
	    jlast = 1;
#line 440 "slatps.f"
	    jinc = -1;
#line 441 "slatps.f"
	}

#line 443 "slatps.f"
	if (tscal != 1.) {
#line 444 "slatps.f"
	    grow = 0.;
#line 445 "slatps.f"
	    goto L80;
#line 446 "slatps.f"
	}

#line 448 "slatps.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 455 "slatps.f"
	    grow = 1. / max(xbnd,smlnum);
#line 456 "slatps.f"
	    xbnd = grow;
#line 457 "slatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 458 "slatps.f"
	    jlen = 1;
#line 459 "slatps.f"
	    i__1 = jlast;
#line 459 "slatps.f"
	    i__2 = jinc;
#line 459 "slatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 463 "slatps.f"
		if (grow <= smlnum) {
#line 463 "slatps.f"
		    goto L80;
#line 463 "slatps.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 468 "slatps.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 469 "slatps.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 469 "slatps.f"
		grow = min(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 473 "slatps.f"
		tjj = (d__1 = ap[ip], abs(d__1));
#line 474 "slatps.f"
		if (xj > tjj) {
#line 474 "slatps.f"
		    xbnd *= tjj / xj;
#line 474 "slatps.f"
		}
#line 476 "slatps.f"
		++jlen;
#line 477 "slatps.f"
		ip += jinc * jlen;
#line 478 "slatps.f"
/* L60: */
#line 478 "slatps.f"
	    }
#line 479 "slatps.f"
	    grow = min(grow,xbnd);
#line 480 "slatps.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 486 "slatps.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 486 "slatps.f"
	    grow = min(d__1,d__2);
#line 487 "slatps.f"
	    i__2 = jlast;
#line 487 "slatps.f"
	    i__1 = jinc;
#line 487 "slatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 491 "slatps.f"
		if (grow <= smlnum) {
#line 491 "slatps.f"
		    goto L80;
#line 491 "slatps.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 496 "slatps.f"
		xj = cnorm[j] + 1.;
#line 497 "slatps.f"
		grow /= xj;
#line 498 "slatps.f"
/* L70: */
#line 498 "slatps.f"
	    }
#line 499 "slatps.f"
	}
#line 500 "slatps.f"
L80:
#line 501 "slatps.f"
	;
#line 501 "slatps.f"
    }

#line 503 "slatps.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 508 "slatps.f"
	stpsv_(uplo, trans, diag, n, &ap[1], &x[1], &c__1, (ftnlen)1, (ftnlen)
		1, (ftnlen)1);
#line 509 "slatps.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 513 "slatps.f"
	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 518 "slatps.f"
	    *scale = bignum / xmax;
#line 519 "slatps.f"
	    sscal_(n, scale, &x[1], &c__1);
#line 520 "slatps.f"
	    xmax = bignum;
#line 521 "slatps.f"
	}

#line 523 "slatps.f"
	if (notran) {

/*           Solve A * x = b */

#line 527 "slatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 528 "slatps.f"
	    i__1 = jlast;
#line 528 "slatps.f"
	    i__2 = jinc;
#line 528 "slatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 532 "slatps.f"
		xj = (d__1 = x[j], abs(d__1));
#line 533 "slatps.f"
		if (nounit) {
#line 534 "slatps.f"
		    tjjs = ap[ip] * tscal;
#line 535 "slatps.f"
		} else {
#line 536 "slatps.f"
		    tjjs = tscal;
#line 537 "slatps.f"
		    if (tscal == 1.) {
#line 537 "slatps.f"
			goto L95;
#line 537 "slatps.f"
		    }
#line 539 "slatps.f"
		}
#line 540 "slatps.f"
		tjj = abs(tjjs);
#line 541 "slatps.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 545 "slatps.f"
		    if (tjj < 1.) {
#line 546 "slatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 550 "slatps.f"
			    rec = 1. / xj;
#line 551 "slatps.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 552 "slatps.f"
			    *scale *= rec;
#line 553 "slatps.f"
			    xmax *= rec;
#line 554 "slatps.f"
			}
#line 555 "slatps.f"
		    }
#line 556 "slatps.f"
		    x[j] /= tjjs;
#line 557 "slatps.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 558 "slatps.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 562 "slatps.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 567 "slatps.f"
			rec = tjj * bignum / xj;
#line 568 "slatps.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 573 "slatps.f"
			    rec /= cnorm[j];
#line 574 "slatps.f"
			}
#line 575 "slatps.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 576 "slatps.f"
			*scale *= rec;
#line 577 "slatps.f"
			xmax *= rec;
#line 578 "slatps.f"
		    }
#line 579 "slatps.f"
		    x[j] /= tjjs;
#line 580 "slatps.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 581 "slatps.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 586 "slatps.f"
		    i__3 = *n;
#line 586 "slatps.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 587 "slatps.f"
			x[i__] = 0.;
#line 588 "slatps.f"
/* L90: */
#line 588 "slatps.f"
		    }
#line 589 "slatps.f"
		    x[j] = 1.;
#line 590 "slatps.f"
		    xj = 1.;
#line 591 "slatps.f"
		    *scale = 0.;
#line 592 "slatps.f"
		    xmax = 0.;
#line 593 "slatps.f"
		}
#line 594 "slatps.f"
L95:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 599 "slatps.f"
		if (xj > 1.) {
#line 600 "slatps.f"
		    rec = 1. / xj;
#line 601 "slatps.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 605 "slatps.f"
			rec *= .5;
#line 606 "slatps.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 607 "slatps.f"
			*scale *= rec;
#line 608 "slatps.f"
		    }
#line 609 "slatps.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 613 "slatps.f"
		    sscal_(n, &c_b36, &x[1], &c__1);
#line 614 "slatps.f"
		    *scale *= .5;
#line 615 "slatps.f"
		}

#line 617 "slatps.f"
		if (upper) {
#line 618 "slatps.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

#line 623 "slatps.f"
			i__3 = j - 1;
#line 623 "slatps.f"
			d__1 = -x[j] * tscal;
#line 623 "slatps.f"
			saxpy_(&i__3, &d__1, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 625 "slatps.f"
			i__3 = j - 1;
#line 625 "slatps.f"
			i__ = isamax_(&i__3, &x[1], &c__1);
#line 626 "slatps.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 627 "slatps.f"
		    }
#line 628 "slatps.f"
		    ip -= j;
#line 629 "slatps.f"
		} else {
#line 630 "slatps.f"
		    if (j < *n) {

/*                    Compute the update */
/*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

#line 635 "slatps.f"
			i__3 = *n - j;
#line 635 "slatps.f"
			d__1 = -x[j] * tscal;
#line 635 "slatps.f"
			saxpy_(&i__3, &d__1, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 637 "slatps.f"
			i__3 = *n - j;
#line 637 "slatps.f"
			i__ = j + isamax_(&i__3, &x[j + 1], &c__1);
#line 638 "slatps.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 639 "slatps.f"
		    }
#line 640 "slatps.f"
		    ip = ip + *n - j + 1;
#line 641 "slatps.f"
		}
#line 642 "slatps.f"
/* L100: */
#line 642 "slatps.f"
	    }

#line 644 "slatps.f"
	} else {

/*           Solve A**T * x = b */

#line 648 "slatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 649 "slatps.f"
	    jlen = 1;
#line 650 "slatps.f"
	    i__2 = jlast;
#line 650 "slatps.f"
	    i__1 = jinc;
#line 650 "slatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 655 "slatps.f"
		xj = (d__1 = x[j], abs(d__1));
#line 656 "slatps.f"
		uscal = tscal;
#line 657 "slatps.f"
		rec = 1. / max(xmax,1.);
#line 658 "slatps.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 662 "slatps.f"
		    rec *= .5;
#line 663 "slatps.f"
		    if (nounit) {
#line 664 "slatps.f"
			tjjs = ap[ip] * tscal;
#line 665 "slatps.f"
		    } else {
#line 666 "slatps.f"
			tjjs = tscal;
#line 667 "slatps.f"
		    }
#line 668 "slatps.f"
		    tjj = abs(tjjs);
#line 669 "slatps.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 673 "slatps.f"
			d__1 = 1., d__2 = rec * tjj;
#line 673 "slatps.f"
			rec = min(d__1,d__2);
#line 674 "slatps.f"
			uscal /= tjjs;
#line 675 "slatps.f"
		    }
#line 676 "slatps.f"
		    if (rec < 1.) {
#line 677 "slatps.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 678 "slatps.f"
			*scale *= rec;
#line 679 "slatps.f"
			xmax *= rec;
#line 680 "slatps.f"
		    }
#line 681 "slatps.f"
		}

#line 683 "slatps.f"
		sumj = 0.;
#line 684 "slatps.f"
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call SDOT to perform the dot product. */

#line 689 "slatps.f"
		    if (upper) {
#line 690 "slatps.f"
			i__3 = j - 1;
#line 690 "slatps.f"
			sumj = sdot_(&i__3, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 691 "slatps.f"
		    } else if (j < *n) {
#line 692 "slatps.f"
			i__3 = *n - j;
#line 692 "slatps.f"
			sumj = sdot_(&i__3, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 693 "slatps.f"
		    }
#line 694 "slatps.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 698 "slatps.f"
		    if (upper) {
#line 699 "slatps.f"
			i__3 = j - 1;
#line 699 "slatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 700 "slatps.f"
			    sumj += ap[ip - j + i__] * uscal * x[i__];
#line 701 "slatps.f"
/* L110: */
#line 701 "slatps.f"
			}
#line 702 "slatps.f"
		    } else if (j < *n) {
#line 703 "slatps.f"
			i__3 = *n - j;
#line 703 "slatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 704 "slatps.f"
			    sumj += ap[ip + i__] * uscal * x[j + i__];
#line 705 "slatps.f"
/* L120: */
#line 705 "slatps.f"
			}
#line 706 "slatps.f"
		    }
#line 707 "slatps.f"
		}

#line 709 "slatps.f"
		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 714 "slatps.f"
		    x[j] -= sumj;
#line 715 "slatps.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 716 "slatps.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 720 "slatps.f"
			tjjs = ap[ip] * tscal;
#line 721 "slatps.f"
		    } else {
#line 722 "slatps.f"
			tjjs = tscal;
#line 723 "slatps.f"
			if (tscal == 1.) {
#line 723 "slatps.f"
			    goto L135;
#line 723 "slatps.f"
			}
#line 725 "slatps.f"
		    }
#line 726 "slatps.f"
		    tjj = abs(tjjs);
#line 727 "slatps.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 731 "slatps.f"
			if (tjj < 1.) {
#line 732 "slatps.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 736 "slatps.f"
				rec = 1. / xj;
#line 737 "slatps.f"
				sscal_(n, &rec, &x[1], &c__1);
#line 738 "slatps.f"
				*scale *= rec;
#line 739 "slatps.f"
				xmax *= rec;
#line 740 "slatps.f"
			    }
#line 741 "slatps.f"
			}
#line 742 "slatps.f"
			x[j] /= tjjs;
#line 743 "slatps.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 747 "slatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 751 "slatps.f"
			    rec = tjj * bignum / xj;
#line 752 "slatps.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 753 "slatps.f"
			    *scale *= rec;
#line 754 "slatps.f"
			    xmax *= rec;
#line 755 "slatps.f"
			}
#line 756 "slatps.f"
			x[j] /= tjjs;
#line 757 "slatps.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0, and compute a solution to A**T*x = 0. */

#line 762 "slatps.f"
			i__3 = *n;
#line 762 "slatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 763 "slatps.f"
			    x[i__] = 0.;
#line 764 "slatps.f"
/* L130: */
#line 764 "slatps.f"
			}
#line 765 "slatps.f"
			x[j] = 1.;
#line 766 "slatps.f"
			*scale = 0.;
#line 767 "slatps.f"
			xmax = 0.;
#line 768 "slatps.f"
		    }
#line 769 "slatps.f"
L135:
#line 770 "slatps.f"
		    ;
#line 770 "slatps.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 775 "slatps.f"
		    x[j] = x[j] / tjjs - sumj;
#line 776 "slatps.f"
		}
/* Computing MAX */
#line 777 "slatps.f"
		d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
#line 777 "slatps.f"
		xmax = max(d__2,d__3);
#line 778 "slatps.f"
		++jlen;
#line 779 "slatps.f"
		ip += jinc * jlen;
#line 780 "slatps.f"
/* L140: */
#line 780 "slatps.f"
	    }
#line 781 "slatps.f"
	}
#line 782 "slatps.f"
	*scale /= tscal;
#line 783 "slatps.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 787 "slatps.f"
    if (tscal != 1.) {
#line 788 "slatps.f"
	d__1 = 1. / tscal;
#line 788 "slatps.f"
	sscal_(n, &d__1, &cnorm[1], &c__1);
#line 789 "slatps.f"
    }

#line 791 "slatps.f"
    return 0;

/*     End of SLATPS */

} /* slatps_ */

