#line 1 "slatrs.f"
/* slatrs.f -- translated by f2c (version 20100827).
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

#line 1 "slatrs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b SLATRS solves a triangular system of equations with the scale factor set to prevent overflow. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLATRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE, */
/*                          CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, LDA, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), CNORM( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLATRS solves one of the triangular systems */
/* > */
/* >    A *x = s*b  or  A**T*x = s*b */
/* > */
/* > with scaling to prevent overflow.  Here A is an upper or lower */
/* > triangular matrix, A**T denotes the transpose of A, x and b are */
/* > n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold.  If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine STRSV is called.  If the matrix A */
/* > is singular (A(j,j) = 0 for some j), then s is set to 0 and a */
/* > non-trivial solution to A*x = 0 is returned. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          The triangular matrix A.  If UPLO = 'U', the leading n by n */
/* >          upper triangular part of the array A contains the upper */
/* >          triangular matrix, and the strictly lower triangular part of */
/* >          A is not referenced.  If UPLO = 'L', the leading n by n lower */
/* >          triangular part of the array A contains the lower triangular */
/* >          matrix, and the strictly upper triangular part of A is not */
/* >          referenced.  If DIAG = 'U', the diagonal elements of A are */
/* >          also not referenced and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max (1,N). */
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

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, STRSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine STRSV if the */
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
/* >  and we can safely call STRSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *a, integer *lda, doublereal *x, 
	doublereal *scale, doublereal *cnorm, integer *info, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len, ftnlen normin_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal xj, rec, tjj;
    static integer jinc;
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
	    integer *, doublereal *, integer *), strsv_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
    static logical notran;
    static integer jfirst;
    static doublereal smlnum;
    static logical nounit;


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 281 "slatrs.f"
    /* Parameter adjustments */
#line 281 "slatrs.f"
    a_dim1 = *lda;
#line 281 "slatrs.f"
    a_offset = 1 + a_dim1;
#line 281 "slatrs.f"
    a -= a_offset;
#line 281 "slatrs.f"
    --x;
#line 281 "slatrs.f"
    --cnorm;
#line 281 "slatrs.f"

#line 281 "slatrs.f"
    /* Function Body */
#line 281 "slatrs.f"
    *info = 0;
#line 282 "slatrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 283 "slatrs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 284 "slatrs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 288 "slatrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 289 "slatrs.f"
	*info = -1;
#line 290 "slatrs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 292 "slatrs.f"
	*info = -2;
#line 293 "slatrs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 294 "slatrs.f"
	*info = -3;
#line 295 "slatrs.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 297 "slatrs.f"
	*info = -4;
#line 298 "slatrs.f"
    } else if (*n < 0) {
#line 299 "slatrs.f"
	*info = -5;
#line 300 "slatrs.f"
    } else if (*lda < max(1,*n)) {
#line 301 "slatrs.f"
	*info = -7;
#line 302 "slatrs.f"
    }
#line 303 "slatrs.f"
    if (*info != 0) {
#line 304 "slatrs.f"
	i__1 = -(*info);
#line 304 "slatrs.f"
	xerbla_("SLATRS", &i__1, (ftnlen)6);
#line 305 "slatrs.f"
	return 0;
#line 306 "slatrs.f"
    }

/*     Quick return if possible */

#line 310 "slatrs.f"
    if (*n == 0) {
#line 310 "slatrs.f"
	return 0;
#line 310 "slatrs.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 315 "slatrs.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 316 "slatrs.f"
    bignum = 1. / smlnum;
#line 317 "slatrs.f"
    *scale = 1.;

#line 319 "slatrs.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 323 "slatrs.f"
	if (upper) {

/*           A is upper triangular. */

#line 327 "slatrs.f"
	    i__1 = *n;
#line 327 "slatrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 328 "slatrs.f"
		i__2 = j - 1;
#line 328 "slatrs.f"
		cnorm[j] = sasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
#line 329 "slatrs.f"
/* L10: */
#line 329 "slatrs.f"
	    }
#line 330 "slatrs.f"
	} else {

/*           A is lower triangular. */

#line 334 "slatrs.f"
	    i__1 = *n - 1;
#line 334 "slatrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 335 "slatrs.f"
		i__2 = *n - j;
#line 335 "slatrs.f"
		cnorm[j] = sasum_(&i__2, &a[j + 1 + j * a_dim1], &c__1);
#line 336 "slatrs.f"
/* L20: */
#line 336 "slatrs.f"
	    }
#line 337 "slatrs.f"
	    cnorm[*n] = 0.;
#line 338 "slatrs.f"
	}
#line 339 "slatrs.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM. */

#line 344 "slatrs.f"
    imax = isamax_(n, &cnorm[1], &c__1);
#line 345 "slatrs.f"
    tmax = cnorm[imax];
#line 346 "slatrs.f"
    if (tmax <= bignum) {
#line 347 "slatrs.f"
	tscal = 1.;
#line 348 "slatrs.f"
    } else {
#line 349 "slatrs.f"
	tscal = 1. / (smlnum * tmax);
#line 350 "slatrs.f"
	sscal_(n, &tscal, &cnorm[1], &c__1);
#line 351 "slatrs.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine STRSV can be used. */

#line 356 "slatrs.f"
    j = isamax_(n, &x[1], &c__1);
#line 357 "slatrs.f"
    xmax = (d__1 = x[j], abs(d__1));
#line 358 "slatrs.f"
    xbnd = xmax;
#line 359 "slatrs.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 363 "slatrs.f"
	if (upper) {
#line 364 "slatrs.f"
	    jfirst = *n;
#line 365 "slatrs.f"
	    jlast = 1;
#line 366 "slatrs.f"
	    jinc = -1;
#line 367 "slatrs.f"
	} else {
#line 368 "slatrs.f"
	    jfirst = 1;
#line 369 "slatrs.f"
	    jlast = *n;
#line 370 "slatrs.f"
	    jinc = 1;
#line 371 "slatrs.f"
	}

#line 373 "slatrs.f"
	if (tscal != 1.) {
#line 374 "slatrs.f"
	    grow = 0.;
#line 375 "slatrs.f"
	    goto L50;
#line 376 "slatrs.f"
	}

#line 378 "slatrs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 385 "slatrs.f"
	    grow = 1. / max(xbnd,smlnum);
#line 386 "slatrs.f"
	    xbnd = grow;
#line 387 "slatrs.f"
	    i__1 = jlast;
#line 387 "slatrs.f"
	    i__2 = jinc;
#line 387 "slatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 391 "slatrs.f"
		if (grow <= smlnum) {
#line 391 "slatrs.f"
		    goto L50;
#line 391 "slatrs.f"
		}

/*              M(j) = G(j-1) / abs(A(j,j)) */

#line 396 "slatrs.f"
		tjj = (d__1 = a[j + j * a_dim1], abs(d__1));
/* Computing MIN */
#line 397 "slatrs.f"
		d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 397 "slatrs.f"
		xbnd = min(d__1,d__2);
#line 398 "slatrs.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 402 "slatrs.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 403 "slatrs.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 407 "slatrs.f"
		    grow = 0.;
#line 408 "slatrs.f"
		}
#line 409 "slatrs.f"
/* L30: */
#line 409 "slatrs.f"
	    }
#line 410 "slatrs.f"
	    grow = xbnd;
#line 411 "slatrs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 417 "slatrs.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 417 "slatrs.f"
	    grow = min(d__1,d__2);
#line 418 "slatrs.f"
	    i__2 = jlast;
#line 418 "slatrs.f"
	    i__1 = jinc;
#line 418 "slatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 422 "slatrs.f"
		if (grow <= smlnum) {
#line 422 "slatrs.f"
		    goto L50;
#line 422 "slatrs.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 427 "slatrs.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 428 "slatrs.f"
/* L40: */
#line 428 "slatrs.f"
	    }
#line 429 "slatrs.f"
	}
#line 430 "slatrs.f"
L50:

#line 432 "slatrs.f"
	;
#line 432 "slatrs.f"
    } else {

/*        Compute the growth in A**T * x = b. */

#line 436 "slatrs.f"
	if (upper) {
#line 437 "slatrs.f"
	    jfirst = 1;
#line 438 "slatrs.f"
	    jlast = *n;
#line 439 "slatrs.f"
	    jinc = 1;
#line 440 "slatrs.f"
	} else {
#line 441 "slatrs.f"
	    jfirst = *n;
#line 442 "slatrs.f"
	    jlast = 1;
#line 443 "slatrs.f"
	    jinc = -1;
#line 444 "slatrs.f"
	}

#line 446 "slatrs.f"
	if (tscal != 1.) {
#line 447 "slatrs.f"
	    grow = 0.;
#line 448 "slatrs.f"
	    goto L80;
#line 449 "slatrs.f"
	}

#line 451 "slatrs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 458 "slatrs.f"
	    grow = 1. / max(xbnd,smlnum);
#line 459 "slatrs.f"
	    xbnd = grow;
#line 460 "slatrs.f"
	    i__1 = jlast;
#line 460 "slatrs.f"
	    i__2 = jinc;
#line 460 "slatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 464 "slatrs.f"
		if (grow <= smlnum) {
#line 464 "slatrs.f"
		    goto L80;
#line 464 "slatrs.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 469 "slatrs.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 470 "slatrs.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 470 "slatrs.f"
		grow = min(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 474 "slatrs.f"
		tjj = (d__1 = a[j + j * a_dim1], abs(d__1));
#line 475 "slatrs.f"
		if (xj > tjj) {
#line 475 "slatrs.f"
		    xbnd *= tjj / xj;
#line 475 "slatrs.f"
		}
#line 477 "slatrs.f"
/* L60: */
#line 477 "slatrs.f"
	    }
#line 478 "slatrs.f"
	    grow = min(grow,xbnd);
#line 479 "slatrs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 485 "slatrs.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 485 "slatrs.f"
	    grow = min(d__1,d__2);
#line 486 "slatrs.f"
	    i__2 = jlast;
#line 486 "slatrs.f"
	    i__1 = jinc;
#line 486 "slatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 490 "slatrs.f"
		if (grow <= smlnum) {
#line 490 "slatrs.f"
		    goto L80;
#line 490 "slatrs.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 495 "slatrs.f"
		xj = cnorm[j] + 1.;
#line 496 "slatrs.f"
		grow /= xj;
#line 497 "slatrs.f"
/* L70: */
#line 497 "slatrs.f"
	    }
#line 498 "slatrs.f"
	}
#line 499 "slatrs.f"
L80:
#line 500 "slatrs.f"
	;
#line 500 "slatrs.f"
    }

#line 502 "slatrs.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 507 "slatrs.f"
	strsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1, (ftnlen)
		1, (ftnlen)1, (ftnlen)1);
#line 508 "slatrs.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 512 "slatrs.f"
	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 517 "slatrs.f"
	    *scale = bignum / xmax;
#line 518 "slatrs.f"
	    sscal_(n, scale, &x[1], &c__1);
#line 519 "slatrs.f"
	    xmax = bignum;
#line 520 "slatrs.f"
	}

#line 522 "slatrs.f"
	if (notran) {

/*           Solve A * x = b */

#line 526 "slatrs.f"
	    i__1 = jlast;
#line 526 "slatrs.f"
	    i__2 = jinc;
#line 526 "slatrs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 530 "slatrs.f"
		xj = (d__1 = x[j], abs(d__1));
#line 531 "slatrs.f"
		if (nounit) {
#line 532 "slatrs.f"
		    tjjs = a[j + j * a_dim1] * tscal;
#line 533 "slatrs.f"
		} else {
#line 534 "slatrs.f"
		    tjjs = tscal;
#line 535 "slatrs.f"
		    if (tscal == 1.) {
#line 535 "slatrs.f"
			goto L95;
#line 535 "slatrs.f"
		    }
#line 537 "slatrs.f"
		}
#line 538 "slatrs.f"
		tjj = abs(tjjs);
#line 539 "slatrs.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 543 "slatrs.f"
		    if (tjj < 1.) {
#line 544 "slatrs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 548 "slatrs.f"
			    rec = 1. / xj;
#line 549 "slatrs.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 550 "slatrs.f"
			    *scale *= rec;
#line 551 "slatrs.f"
			    xmax *= rec;
#line 552 "slatrs.f"
			}
#line 553 "slatrs.f"
		    }
#line 554 "slatrs.f"
		    x[j] /= tjjs;
#line 555 "slatrs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 556 "slatrs.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 560 "slatrs.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 565 "slatrs.f"
			rec = tjj * bignum / xj;
#line 566 "slatrs.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 571 "slatrs.f"
			    rec /= cnorm[j];
#line 572 "slatrs.f"
			}
#line 573 "slatrs.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 574 "slatrs.f"
			*scale *= rec;
#line 575 "slatrs.f"
			xmax *= rec;
#line 576 "slatrs.f"
		    }
#line 577 "slatrs.f"
		    x[j] /= tjjs;
#line 578 "slatrs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 579 "slatrs.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 584 "slatrs.f"
		    i__3 = *n;
#line 584 "slatrs.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 585 "slatrs.f"
			x[i__] = 0.;
#line 586 "slatrs.f"
/* L90: */
#line 586 "slatrs.f"
		    }
#line 587 "slatrs.f"
		    x[j] = 1.;
#line 588 "slatrs.f"
		    xj = 1.;
#line 589 "slatrs.f"
		    *scale = 0.;
#line 590 "slatrs.f"
		    xmax = 0.;
#line 591 "slatrs.f"
		}
#line 592 "slatrs.f"
L95:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 597 "slatrs.f"
		if (xj > 1.) {
#line 598 "slatrs.f"
		    rec = 1. / xj;
#line 599 "slatrs.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 603 "slatrs.f"
			rec *= .5;
#line 604 "slatrs.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 605 "slatrs.f"
			*scale *= rec;
#line 606 "slatrs.f"
		    }
#line 607 "slatrs.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 611 "slatrs.f"
		    sscal_(n, &c_b36, &x[1], &c__1);
#line 612 "slatrs.f"
		    *scale *= .5;
#line 613 "slatrs.f"
		}

#line 615 "slatrs.f"
		if (upper) {
#line 616 "slatrs.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

#line 621 "slatrs.f"
			i__3 = j - 1;
#line 621 "slatrs.f"
			d__1 = -x[j] * tscal;
#line 621 "slatrs.f"
			saxpy_(&i__3, &d__1, &a[j * a_dim1 + 1], &c__1, &x[1],
				 &c__1);
#line 623 "slatrs.f"
			i__3 = j - 1;
#line 623 "slatrs.f"
			i__ = isamax_(&i__3, &x[1], &c__1);
#line 624 "slatrs.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 625 "slatrs.f"
		    }
#line 626 "slatrs.f"
		} else {
#line 627 "slatrs.f"
		    if (j < *n) {

/*                    Compute the update */
/*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

#line 632 "slatrs.f"
			i__3 = *n - j;
#line 632 "slatrs.f"
			d__1 = -x[j] * tscal;
#line 632 "slatrs.f"
			saxpy_(&i__3, &d__1, &a[j + 1 + j * a_dim1], &c__1, &
				x[j + 1], &c__1);
#line 634 "slatrs.f"
			i__3 = *n - j;
#line 634 "slatrs.f"
			i__ = j + isamax_(&i__3, &x[j + 1], &c__1);
#line 635 "slatrs.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 636 "slatrs.f"
		    }
#line 637 "slatrs.f"
		}
#line 638 "slatrs.f"
/* L100: */
#line 638 "slatrs.f"
	    }

#line 640 "slatrs.f"
	} else {

/*           Solve A**T * x = b */

#line 644 "slatrs.f"
	    i__2 = jlast;
#line 644 "slatrs.f"
	    i__1 = jinc;
#line 644 "slatrs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 649 "slatrs.f"
		xj = (d__1 = x[j], abs(d__1));
#line 650 "slatrs.f"
		uscal = tscal;
#line 651 "slatrs.f"
		rec = 1. / max(xmax,1.);
#line 652 "slatrs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 656 "slatrs.f"
		    rec *= .5;
#line 657 "slatrs.f"
		    if (nounit) {
#line 658 "slatrs.f"
			tjjs = a[j + j * a_dim1] * tscal;
#line 659 "slatrs.f"
		    } else {
#line 660 "slatrs.f"
			tjjs = tscal;
#line 661 "slatrs.f"
		    }
#line 662 "slatrs.f"
		    tjj = abs(tjjs);
#line 663 "slatrs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 667 "slatrs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 667 "slatrs.f"
			rec = min(d__1,d__2);
#line 668 "slatrs.f"
			uscal /= tjjs;
#line 669 "slatrs.f"
		    }
#line 670 "slatrs.f"
		    if (rec < 1.) {
#line 671 "slatrs.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 672 "slatrs.f"
			*scale *= rec;
#line 673 "slatrs.f"
			xmax *= rec;
#line 674 "slatrs.f"
		    }
#line 675 "slatrs.f"
		}

#line 677 "slatrs.f"
		sumj = 0.;
#line 678 "slatrs.f"
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call SDOT to perform the dot product. */

#line 683 "slatrs.f"
		    if (upper) {
#line 684 "slatrs.f"
			i__3 = j - 1;
#line 684 "slatrs.f"
			sumj = sdot_(&i__3, &a[j * a_dim1 + 1], &c__1, &x[1], 
				&c__1);
#line 685 "slatrs.f"
		    } else if (j < *n) {
#line 686 "slatrs.f"
			i__3 = *n - j;
#line 686 "slatrs.f"
			sumj = sdot_(&i__3, &a[j + 1 + j * a_dim1], &c__1, &x[
				j + 1], &c__1);
#line 687 "slatrs.f"
		    }
#line 688 "slatrs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 692 "slatrs.f"
		    if (upper) {
#line 693 "slatrs.f"
			i__3 = j - 1;
#line 693 "slatrs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 694 "slatrs.f"
			    sumj += a[i__ + j * a_dim1] * uscal * x[i__];
#line 695 "slatrs.f"
/* L110: */
#line 695 "slatrs.f"
			}
#line 696 "slatrs.f"
		    } else if (j < *n) {
#line 697 "slatrs.f"
			i__3 = *n;
#line 697 "slatrs.f"
			for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 698 "slatrs.f"
			    sumj += a[i__ + j * a_dim1] * uscal * x[i__];
#line 699 "slatrs.f"
/* L120: */
#line 699 "slatrs.f"
			}
#line 700 "slatrs.f"
		    }
#line 701 "slatrs.f"
		}

#line 703 "slatrs.f"
		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 708 "slatrs.f"
		    x[j] -= sumj;
#line 709 "slatrs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 710 "slatrs.f"
		    if (nounit) {
#line 711 "slatrs.f"
			tjjs = a[j + j * a_dim1] * tscal;
#line 712 "slatrs.f"
		    } else {
#line 713 "slatrs.f"
			tjjs = tscal;
#line 714 "slatrs.f"
			if (tscal == 1.) {
#line 714 "slatrs.f"
			    goto L135;
#line 714 "slatrs.f"
			}
#line 716 "slatrs.f"
		    }

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 720 "slatrs.f"
		    tjj = abs(tjjs);
#line 721 "slatrs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 725 "slatrs.f"
			if (tjj < 1.) {
#line 726 "slatrs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 730 "slatrs.f"
				rec = 1. / xj;
#line 731 "slatrs.f"
				sscal_(n, &rec, &x[1], &c__1);
#line 732 "slatrs.f"
				*scale *= rec;
#line 733 "slatrs.f"
				xmax *= rec;
#line 734 "slatrs.f"
			    }
#line 735 "slatrs.f"
			}
#line 736 "slatrs.f"
			x[j] /= tjjs;
#line 737 "slatrs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 741 "slatrs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 745 "slatrs.f"
			    rec = tjj * bignum / xj;
#line 746 "slatrs.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 747 "slatrs.f"
			    *scale *= rec;
#line 748 "slatrs.f"
			    xmax *= rec;
#line 749 "slatrs.f"
			}
#line 750 "slatrs.f"
			x[j] /= tjjs;
#line 751 "slatrs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0, and compute a solution to A**T*x = 0. */

#line 756 "slatrs.f"
			i__3 = *n;
#line 756 "slatrs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 757 "slatrs.f"
			    x[i__] = 0.;
#line 758 "slatrs.f"
/* L130: */
#line 758 "slatrs.f"
			}
#line 759 "slatrs.f"
			x[j] = 1.;
#line 760 "slatrs.f"
			*scale = 0.;
#line 761 "slatrs.f"
			xmax = 0.;
#line 762 "slatrs.f"
		    }
#line 763 "slatrs.f"
L135:
#line 764 "slatrs.f"
		    ;
#line 764 "slatrs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 769 "slatrs.f"
		    x[j] = x[j] / tjjs - sumj;
#line 770 "slatrs.f"
		}
/* Computing MAX */
#line 771 "slatrs.f"
		d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
#line 771 "slatrs.f"
		xmax = max(d__2,d__3);
#line 772 "slatrs.f"
/* L140: */
#line 772 "slatrs.f"
	    }
#line 773 "slatrs.f"
	}
#line 774 "slatrs.f"
	*scale /= tscal;
#line 775 "slatrs.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 779 "slatrs.f"
    if (tscal != 1.) {
#line 780 "slatrs.f"
	d__1 = 1. / tscal;
#line 780 "slatrs.f"
	sscal_(n, &d__1, &cnorm[1], &c__1);
#line 781 "slatrs.f"
    }

#line 783 "slatrs.f"
    return 0;

/*     End of SLATRS */

} /* slatrs_ */

