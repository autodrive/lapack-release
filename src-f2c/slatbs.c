#line 1 "slatbs.f"
/* slatbs.f -- translated by f2c (version 20100827).
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

#line 1 "slatbs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b SLATBS solves a triangular banded system of equations. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLATBS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slatbs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slatbs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slatbs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, */
/*                          SCALE, CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), CNORM( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLATBS solves one of the triangular systems */
/* > */
/* >    A *x = s*b  or  A**T*x = s*b */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular band matrix.  Here A**T denotes the transpose of A, x and b */
/* > are n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold.  If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine STBSV is called.  If the matrix A */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of subdiagonals or superdiagonals in the */
/* >          triangular matrix A.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          The upper or lower triangular band matrix A, stored in the */
/* >          first KD+1 rows of the array. The j-th column of A is stored */
/* >          in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
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
/* >  A rough bound on x is computed; if that is less than overflow, STBSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine STBSV if the */
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
/* >  and we can safely call STBSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublereal *ab, integer *ldab, 
	doublereal *x, doublereal *scale, doublereal *cnorm, integer *info, 
	ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len, ftnlen normin_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal xj, rec, tjj;
    static integer jinc, jlen;
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax, tjjs;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal xmax, grow, sumj;
    static integer maind;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal tscal, uscal;
    static integer jlast;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int stbsv_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
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

#line 285 "slatbs.f"
    /* Parameter adjustments */
#line 285 "slatbs.f"
    ab_dim1 = *ldab;
#line 285 "slatbs.f"
    ab_offset = 1 + ab_dim1;
#line 285 "slatbs.f"
    ab -= ab_offset;
#line 285 "slatbs.f"
    --x;
#line 285 "slatbs.f"
    --cnorm;
#line 285 "slatbs.f"

#line 285 "slatbs.f"
    /* Function Body */
#line 285 "slatbs.f"
    *info = 0;
#line 286 "slatbs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 287 "slatbs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 288 "slatbs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 292 "slatbs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 293 "slatbs.f"
	*info = -1;
#line 294 "slatbs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 296 "slatbs.f"
	*info = -2;
#line 297 "slatbs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 298 "slatbs.f"
	*info = -3;
#line 299 "slatbs.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 301 "slatbs.f"
	*info = -4;
#line 302 "slatbs.f"
    } else if (*n < 0) {
#line 303 "slatbs.f"
	*info = -5;
#line 304 "slatbs.f"
    } else if (*kd < 0) {
#line 305 "slatbs.f"
	*info = -6;
#line 306 "slatbs.f"
    } else if (*ldab < *kd + 1) {
#line 307 "slatbs.f"
	*info = -8;
#line 308 "slatbs.f"
    }
#line 309 "slatbs.f"
    if (*info != 0) {
#line 310 "slatbs.f"
	i__1 = -(*info);
#line 310 "slatbs.f"
	xerbla_("SLATBS", &i__1, (ftnlen)6);
#line 311 "slatbs.f"
	return 0;
#line 312 "slatbs.f"
    }

/*     Quick return if possible */

#line 316 "slatbs.f"
    if (*n == 0) {
#line 316 "slatbs.f"
	return 0;
#line 316 "slatbs.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 321 "slatbs.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12) / slamch_("Precision", (
	    ftnlen)9);
#line 322 "slatbs.f"
    bignum = 1. / smlnum;
#line 323 "slatbs.f"
    *scale = 1.;

#line 325 "slatbs.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 329 "slatbs.f"
	if (upper) {

/*           A is upper triangular. */

#line 333 "slatbs.f"
	    i__1 = *n;
#line 333 "slatbs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 334 "slatbs.f"
		i__2 = *kd, i__3 = j - 1;
#line 334 "slatbs.f"
		jlen = min(i__2,i__3);
#line 335 "slatbs.f"
		cnorm[j] = sasum_(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1], &
			c__1);
#line 336 "slatbs.f"
/* L10: */
#line 336 "slatbs.f"
	    }
#line 337 "slatbs.f"
	} else {

/*           A is lower triangular. */

#line 341 "slatbs.f"
	    i__1 = *n;
#line 341 "slatbs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 342 "slatbs.f"
		i__2 = *kd, i__3 = *n - j;
#line 342 "slatbs.f"
		jlen = min(i__2,i__3);
#line 343 "slatbs.f"
		if (jlen > 0) {
#line 344 "slatbs.f"
		    cnorm[j] = sasum_(&jlen, &ab[j * ab_dim1 + 2], &c__1);
#line 345 "slatbs.f"
		} else {
#line 346 "slatbs.f"
		    cnorm[j] = 0.;
#line 347 "slatbs.f"
		}
#line 348 "slatbs.f"
/* L20: */
#line 348 "slatbs.f"
	    }
#line 349 "slatbs.f"
	}
#line 350 "slatbs.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM. */

#line 355 "slatbs.f"
    imax = isamax_(n, &cnorm[1], &c__1);
#line 356 "slatbs.f"
    tmax = cnorm[imax];
#line 357 "slatbs.f"
    if (tmax <= bignum) {
#line 358 "slatbs.f"
	tscal = 1.;
#line 359 "slatbs.f"
    } else {
#line 360 "slatbs.f"
	tscal = 1. / (smlnum * tmax);
#line 361 "slatbs.f"
	sscal_(n, &tscal, &cnorm[1], &c__1);
#line 362 "slatbs.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine STBSV can be used. */

#line 367 "slatbs.f"
    j = isamax_(n, &x[1], &c__1);
#line 368 "slatbs.f"
    xmax = (d__1 = x[j], abs(d__1));
#line 369 "slatbs.f"
    xbnd = xmax;
#line 370 "slatbs.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 374 "slatbs.f"
	if (upper) {
#line 375 "slatbs.f"
	    jfirst = *n;
#line 376 "slatbs.f"
	    jlast = 1;
#line 377 "slatbs.f"
	    jinc = -1;
#line 378 "slatbs.f"
	    maind = *kd + 1;
#line 379 "slatbs.f"
	} else {
#line 380 "slatbs.f"
	    jfirst = 1;
#line 381 "slatbs.f"
	    jlast = *n;
#line 382 "slatbs.f"
	    jinc = 1;
#line 383 "slatbs.f"
	    maind = 1;
#line 384 "slatbs.f"
	}

#line 386 "slatbs.f"
	if (tscal != 1.) {
#line 387 "slatbs.f"
	    grow = 0.;
#line 388 "slatbs.f"
	    goto L50;
#line 389 "slatbs.f"
	}

#line 391 "slatbs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 398 "slatbs.f"
	    grow = 1. / max(xbnd,smlnum);
#line 399 "slatbs.f"
	    xbnd = grow;
#line 400 "slatbs.f"
	    i__1 = jlast;
#line 400 "slatbs.f"
	    i__2 = jinc;
#line 400 "slatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 404 "slatbs.f"
		if (grow <= smlnum) {
#line 404 "slatbs.f"
		    goto L50;
#line 404 "slatbs.f"
		}

/*              M(j) = G(j-1) / abs(A(j,j)) */

#line 409 "slatbs.f"
		tjj = (d__1 = ab[maind + j * ab_dim1], abs(d__1));
/* Computing MIN */
#line 410 "slatbs.f"
		d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 410 "slatbs.f"
		xbnd = min(d__1,d__2);
#line 411 "slatbs.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 415 "slatbs.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 416 "slatbs.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 420 "slatbs.f"
		    grow = 0.;
#line 421 "slatbs.f"
		}
#line 422 "slatbs.f"
/* L30: */
#line 422 "slatbs.f"
	    }
#line 423 "slatbs.f"
	    grow = xbnd;
#line 424 "slatbs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 430 "slatbs.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 430 "slatbs.f"
	    grow = min(d__1,d__2);
#line 431 "slatbs.f"
	    i__2 = jlast;
#line 431 "slatbs.f"
	    i__1 = jinc;
#line 431 "slatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 435 "slatbs.f"
		if (grow <= smlnum) {
#line 435 "slatbs.f"
		    goto L50;
#line 435 "slatbs.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 440 "slatbs.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 441 "slatbs.f"
/* L40: */
#line 441 "slatbs.f"
	    }
#line 442 "slatbs.f"
	}
#line 443 "slatbs.f"
L50:

#line 445 "slatbs.f"
	;
#line 445 "slatbs.f"
    } else {

/*        Compute the growth in A**T * x = b. */

#line 449 "slatbs.f"
	if (upper) {
#line 450 "slatbs.f"
	    jfirst = 1;
#line 451 "slatbs.f"
	    jlast = *n;
#line 452 "slatbs.f"
	    jinc = 1;
#line 453 "slatbs.f"
	    maind = *kd + 1;
#line 454 "slatbs.f"
	} else {
#line 455 "slatbs.f"
	    jfirst = *n;
#line 456 "slatbs.f"
	    jlast = 1;
#line 457 "slatbs.f"
	    jinc = -1;
#line 458 "slatbs.f"
	    maind = 1;
#line 459 "slatbs.f"
	}

#line 461 "slatbs.f"
	if (tscal != 1.) {
#line 462 "slatbs.f"
	    grow = 0.;
#line 463 "slatbs.f"
	    goto L80;
#line 464 "slatbs.f"
	}

#line 466 "slatbs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 473 "slatbs.f"
	    grow = 1. / max(xbnd,smlnum);
#line 474 "slatbs.f"
	    xbnd = grow;
#line 475 "slatbs.f"
	    i__1 = jlast;
#line 475 "slatbs.f"
	    i__2 = jinc;
#line 475 "slatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 479 "slatbs.f"
		if (grow <= smlnum) {
#line 479 "slatbs.f"
		    goto L80;
#line 479 "slatbs.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 484 "slatbs.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 485 "slatbs.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 485 "slatbs.f"
		grow = min(d__1,d__2);

/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 489 "slatbs.f"
		tjj = (d__1 = ab[maind + j * ab_dim1], abs(d__1));
#line 490 "slatbs.f"
		if (xj > tjj) {
#line 490 "slatbs.f"
		    xbnd *= tjj / xj;
#line 490 "slatbs.f"
		}
#line 492 "slatbs.f"
/* L60: */
#line 492 "slatbs.f"
	    }
#line 493 "slatbs.f"
	    grow = min(grow,xbnd);
#line 494 "slatbs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 500 "slatbs.f"
	    d__1 = 1., d__2 = 1. / max(xbnd,smlnum);
#line 500 "slatbs.f"
	    grow = min(d__1,d__2);
#line 501 "slatbs.f"
	    i__2 = jlast;
#line 501 "slatbs.f"
	    i__1 = jinc;
#line 501 "slatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 505 "slatbs.f"
		if (grow <= smlnum) {
#line 505 "slatbs.f"
		    goto L80;
#line 505 "slatbs.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 510 "slatbs.f"
		xj = cnorm[j] + 1.;
#line 511 "slatbs.f"
		grow /= xj;
#line 512 "slatbs.f"
/* L70: */
#line 512 "slatbs.f"
	    }
#line 513 "slatbs.f"
	}
#line 514 "slatbs.f"
L80:
#line 515 "slatbs.f"
	;
#line 515 "slatbs.f"
    }

#line 517 "slatbs.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 522 "slatbs.f"
	stbsv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &x[1], &c__1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 523 "slatbs.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 527 "slatbs.f"
	if (xmax > bignum) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 532 "slatbs.f"
	    *scale = bignum / xmax;
#line 533 "slatbs.f"
	    sscal_(n, scale, &x[1], &c__1);
#line 534 "slatbs.f"
	    xmax = bignum;
#line 535 "slatbs.f"
	}

#line 537 "slatbs.f"
	if (notran) {

/*           Solve A * x = b */

#line 541 "slatbs.f"
	    i__1 = jlast;
#line 541 "slatbs.f"
	    i__2 = jinc;
#line 541 "slatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 545 "slatbs.f"
		xj = (d__1 = x[j], abs(d__1));
#line 546 "slatbs.f"
		if (nounit) {
#line 547 "slatbs.f"
		    tjjs = ab[maind + j * ab_dim1] * tscal;
#line 548 "slatbs.f"
		} else {
#line 549 "slatbs.f"
		    tjjs = tscal;
#line 550 "slatbs.f"
		    if (tscal == 1.) {
#line 550 "slatbs.f"
			goto L95;
#line 550 "slatbs.f"
		    }
#line 552 "slatbs.f"
		}
#line 553 "slatbs.f"
		tjj = abs(tjjs);
#line 554 "slatbs.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 558 "slatbs.f"
		    if (tjj < 1.) {
#line 559 "slatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 563 "slatbs.f"
			    rec = 1. / xj;
#line 564 "slatbs.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 565 "slatbs.f"
			    *scale *= rec;
#line 566 "slatbs.f"
			    xmax *= rec;
#line 567 "slatbs.f"
			}
#line 568 "slatbs.f"
		    }
#line 569 "slatbs.f"
		    x[j] /= tjjs;
#line 570 "slatbs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 571 "slatbs.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 575 "slatbs.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 580 "slatbs.f"
			rec = tjj * bignum / xj;
#line 581 "slatbs.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 586 "slatbs.f"
			    rec /= cnorm[j];
#line 587 "slatbs.f"
			}
#line 588 "slatbs.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 589 "slatbs.f"
			*scale *= rec;
#line 590 "slatbs.f"
			xmax *= rec;
#line 591 "slatbs.f"
		    }
#line 592 "slatbs.f"
		    x[j] /= tjjs;
#line 593 "slatbs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 594 "slatbs.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 599 "slatbs.f"
		    i__3 = *n;
#line 599 "slatbs.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 600 "slatbs.f"
			x[i__] = 0.;
#line 601 "slatbs.f"
/* L90: */
#line 601 "slatbs.f"
		    }
#line 602 "slatbs.f"
		    x[j] = 1.;
#line 603 "slatbs.f"
		    xj = 1.;
#line 604 "slatbs.f"
		    *scale = 0.;
#line 605 "slatbs.f"
		    xmax = 0.;
#line 606 "slatbs.f"
		}
#line 607 "slatbs.f"
L95:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 612 "slatbs.f"
		if (xj > 1.) {
#line 613 "slatbs.f"
		    rec = 1. / xj;
#line 614 "slatbs.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 618 "slatbs.f"
			rec *= .5;
#line 619 "slatbs.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 620 "slatbs.f"
			*scale *= rec;
#line 621 "slatbs.f"
		    }
#line 622 "slatbs.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 626 "slatbs.f"
		    sscal_(n, &c_b36, &x[1], &c__1);
#line 627 "slatbs.f"
		    *scale *= .5;
#line 628 "slatbs.f"
		}

#line 630 "slatbs.f"
		if (upper) {
#line 631 "slatbs.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) - */
/*                                             x(j)* A(max(1,j-kd):j-1,j) */

/* Computing MIN */
#line 637 "slatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 637 "slatbs.f"
			jlen = min(i__3,i__4);
#line 638 "slatbs.f"
			d__1 = -x[j] * tscal;
#line 638 "slatbs.f"
			saxpy_(&jlen, &d__1, &ab[*kd + 1 - jlen + j * ab_dim1]
				, &c__1, &x[j - jlen], &c__1);
#line 640 "slatbs.f"
			i__3 = j - 1;
#line 640 "slatbs.f"
			i__ = isamax_(&i__3, &x[1], &c__1);
#line 641 "slatbs.f"
			xmax = (d__1 = x[i__], abs(d__1));
#line 642 "slatbs.f"
		    }
#line 643 "slatbs.f"
		} else if (j < *n) {

/*                 Compute the update */
/*                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) - */
/*                                          x(j) * A(j+1:min(j+kd,n),j) */

/* Computing MIN */
#line 649 "slatbs.f"
		    i__3 = *kd, i__4 = *n - j;
#line 649 "slatbs.f"
		    jlen = min(i__3,i__4);
#line 650 "slatbs.f"
		    if (jlen > 0) {
#line 650 "slatbs.f"
			d__1 = -x[j] * tscal;
#line 650 "slatbs.f"
			saxpy_(&jlen, &d__1, &ab[j * ab_dim1 + 2], &c__1, &x[
				j + 1], &c__1);
#line 650 "slatbs.f"
		    }
#line 653 "slatbs.f"
		    i__3 = *n - j;
#line 653 "slatbs.f"
		    i__ = j + isamax_(&i__3, &x[j + 1], &c__1);
#line 654 "slatbs.f"
		    xmax = (d__1 = x[i__], abs(d__1));
#line 655 "slatbs.f"
		}
#line 656 "slatbs.f"
/* L100: */
#line 656 "slatbs.f"
	    }

#line 658 "slatbs.f"
	} else {

/*           Solve A**T * x = b */

#line 662 "slatbs.f"
	    i__2 = jlast;
#line 662 "slatbs.f"
	    i__1 = jinc;
#line 662 "slatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 667 "slatbs.f"
		xj = (d__1 = x[j], abs(d__1));
#line 668 "slatbs.f"
		uscal = tscal;
#line 669 "slatbs.f"
		rec = 1. / max(xmax,1.);
#line 670 "slatbs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 674 "slatbs.f"
		    rec *= .5;
#line 675 "slatbs.f"
		    if (nounit) {
#line 676 "slatbs.f"
			tjjs = ab[maind + j * ab_dim1] * tscal;
#line 677 "slatbs.f"
		    } else {
#line 678 "slatbs.f"
			tjjs = tscal;
#line 679 "slatbs.f"
		    }
#line 680 "slatbs.f"
		    tjj = abs(tjjs);
#line 681 "slatbs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 685 "slatbs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 685 "slatbs.f"
			rec = min(d__1,d__2);
#line 686 "slatbs.f"
			uscal /= tjjs;
#line 687 "slatbs.f"
		    }
#line 688 "slatbs.f"
		    if (rec < 1.) {
#line 689 "slatbs.f"
			sscal_(n, &rec, &x[1], &c__1);
#line 690 "slatbs.f"
			*scale *= rec;
#line 691 "slatbs.f"
			xmax *= rec;
#line 692 "slatbs.f"
		    }
#line 693 "slatbs.f"
		}

#line 695 "slatbs.f"
		sumj = 0.;
#line 696 "slatbs.f"
		if (uscal == 1.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call SDOT to perform the dot product. */

#line 701 "slatbs.f"
		    if (upper) {
/* Computing MIN */
#line 702 "slatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 702 "slatbs.f"
			jlen = min(i__3,i__4);
#line 703 "slatbs.f"
			sumj = sdot_(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1],
				 &c__1, &x[j - jlen], &c__1);
#line 705 "slatbs.f"
		    } else {
/* Computing MIN */
#line 706 "slatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 706 "slatbs.f"
			jlen = min(i__3,i__4);
#line 707 "slatbs.f"
			if (jlen > 0) {
#line 707 "slatbs.f"
			    sumj = sdot_(&jlen, &ab[j * ab_dim1 + 2], &c__1, &
				    x[j + 1], &c__1);
#line 707 "slatbs.f"
			}
#line 709 "slatbs.f"
		    }
#line 710 "slatbs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 714 "slatbs.f"
		    if (upper) {
/* Computing MIN */
#line 715 "slatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 715 "slatbs.f"
			jlen = min(i__3,i__4);
#line 716 "slatbs.f"
			i__3 = jlen;
#line 716 "slatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 717 "slatbs.f"
			    sumj += ab[*kd + i__ - jlen + j * ab_dim1] * 
				    uscal * x[j - jlen - 1 + i__];
#line 719 "slatbs.f"
/* L110: */
#line 719 "slatbs.f"
			}
#line 720 "slatbs.f"
		    } else {
/* Computing MIN */
#line 721 "slatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 721 "slatbs.f"
			jlen = min(i__3,i__4);
#line 722 "slatbs.f"
			i__3 = jlen;
#line 722 "slatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 723 "slatbs.f"
			    sumj += ab[i__ + 1 + j * ab_dim1] * uscal * x[j + 
				    i__];
#line 724 "slatbs.f"
/* L120: */
#line 724 "slatbs.f"
			}
#line 725 "slatbs.f"
		    }
#line 726 "slatbs.f"
		}

#line 728 "slatbs.f"
		if (uscal == tscal) {

/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 733 "slatbs.f"
		    x[j] -= sumj;
#line 734 "slatbs.f"
		    xj = (d__1 = x[j], abs(d__1));
#line 735 "slatbs.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 739 "slatbs.f"
			tjjs = ab[maind + j * ab_dim1] * tscal;
#line 740 "slatbs.f"
		    } else {
#line 741 "slatbs.f"
			tjjs = tscal;
#line 742 "slatbs.f"
			if (tscal == 1.) {
#line 742 "slatbs.f"
			    goto L135;
#line 742 "slatbs.f"
			}
#line 744 "slatbs.f"
		    }
#line 745 "slatbs.f"
		    tjj = abs(tjjs);
#line 746 "slatbs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 750 "slatbs.f"
			if (tjj < 1.) {
#line 751 "slatbs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 755 "slatbs.f"
				rec = 1. / xj;
#line 756 "slatbs.f"
				sscal_(n, &rec, &x[1], &c__1);
#line 757 "slatbs.f"
				*scale *= rec;
#line 758 "slatbs.f"
				xmax *= rec;
#line 759 "slatbs.f"
			    }
#line 760 "slatbs.f"
			}
#line 761 "slatbs.f"
			x[j] /= tjjs;
#line 762 "slatbs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 766 "slatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 770 "slatbs.f"
			    rec = tjj * bignum / xj;
#line 771 "slatbs.f"
			    sscal_(n, &rec, &x[1], &c__1);
#line 772 "slatbs.f"
			    *scale *= rec;
#line 773 "slatbs.f"
			    xmax *= rec;
#line 774 "slatbs.f"
			}
#line 775 "slatbs.f"
			x[j] /= tjjs;
#line 776 "slatbs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0, and compute a solution to A**T*x = 0. */

#line 781 "slatbs.f"
			i__3 = *n;
#line 781 "slatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 782 "slatbs.f"
			    x[i__] = 0.;
#line 783 "slatbs.f"
/* L130: */
#line 783 "slatbs.f"
			}
#line 784 "slatbs.f"
			x[j] = 1.;
#line 785 "slatbs.f"
			*scale = 0.;
#line 786 "slatbs.f"
			xmax = 0.;
#line 787 "slatbs.f"
		    }
#line 788 "slatbs.f"
L135:
#line 789 "slatbs.f"
		    ;
#line 789 "slatbs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - sumj if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 794 "slatbs.f"
		    x[j] = x[j] / tjjs - sumj;
#line 795 "slatbs.f"
		}
/* Computing MAX */
#line 796 "slatbs.f"
		d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
#line 796 "slatbs.f"
		xmax = max(d__2,d__3);
#line 797 "slatbs.f"
/* L140: */
#line 797 "slatbs.f"
	    }
#line 798 "slatbs.f"
	}
#line 799 "slatbs.f"
	*scale /= tscal;
#line 800 "slatbs.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 804 "slatbs.f"
    if (tscal != 1.) {
#line 805 "slatbs.f"
	d__1 = 1. / tscal;
#line 805 "slatbs.f"
	sscal_(n, &d__1, &cnorm[1], &c__1);
#line 806 "slatbs.f"
    }

#line 808 "slatbs.f"
    return 0;

/*     End of SLATBS */

} /* slatbs_ */

