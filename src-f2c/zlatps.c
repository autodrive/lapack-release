#line 1 "zlatps.f"
/* zlatps.f -- translated by f2c (version 20100827).
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

#line 1 "zlatps.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b ZLATPS solves a triangular system of equations with the matrix held in packed storage. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLATPS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatps.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatps.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatps.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLATPS( UPLO, TRANS, DIAG, NORMIN, N, AP, X, SCALE, */
/*                          CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   CNORM( * ) */
/*       COMPLEX*16         AP( * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLATPS solves one of the triangular systems */
/* > */
/* >    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b, */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular matrix stored in packed form.  Here A**T denotes the */
/* > transpose of A, A**H denotes the conjugate transpose of A, x and b */
/* > are n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold.  If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine ZTPSV is called. If the matrix A */
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
/* >          = 'N':  Solve A * x = s*b     (No transpose) */
/* >          = 'T':  Solve A**T * x = s*b  (Transpose) */
/* >          = 'C':  Solve A**H * x = s*b  (Conjugate transpose) */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (N) */
/* >          On entry, the right hand side b of the triangular system. */
/* >          On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          The scaling factor s for the triangular system */
/* >             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b. */
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

/* > \ingroup complex16OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, ZTPSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTPSV if the */
/* >  reciprocal of the largest M(j), j=1,..,n, is larger than */
/* >  max(underflow, 1/overflow). */
/* > */
/* >  The bound on x(j) is also used to determine when a step in the */
/* >  columnwise method can be performed without fear of overflow.  If */
/* >  the computed bound is greater than a large constant, x is scaled to */
/* >  prevent overflow, but if the bound overflows, x is set to 0, x(j) to */
/* >  1, and scale to 0, and a non-trivial solution to A*x = 0 is found. */
/* > */
/* >  Similarly, a row-wise scheme is used to solve A**T *x = b  or */
/* >  A**H *x = b.  The basic algorithm for A upper triangular is */
/* > */
/* >       for j = 1, ..., n */
/* >            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j) */
/* >       end */
/* > */
/* >  We simultaneously compute two bounds */
/* >       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j */
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
/* >  and we can safely call ZTPSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *ap, doublecomplex *x, doublereal *
	scale, doublereal *cnorm, integer *info, ftnlen uplo_len, ftnlen 
	trans_len, ftnlen diag_len, ftnlen normin_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, ip;
    static doublereal xj, rec, tjj;
    static integer jinc, jlen;
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax;
    static doublecomplex tjjs;
    static doublereal xmax, grow;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal tscal;
    static doublecomplex uscal;
    static integer jlast;
    static doublecomplex csumj;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Double Complex */ VOID zdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ztpsv_(
	    char *, char *, char *, integer *, doublecomplex *, doublecomplex 
	    *, integer *, ftnlen, ftnlen, ftnlen), dlabad_(doublereal *, 
	    doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal bignum;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static logical notran;
    static integer jfirst;
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 287 "zlatps.f"
    /* Parameter adjustments */
#line 287 "zlatps.f"
    --cnorm;
#line 287 "zlatps.f"
    --x;
#line 287 "zlatps.f"
    --ap;
#line 287 "zlatps.f"

#line 287 "zlatps.f"
    /* Function Body */
#line 287 "zlatps.f"
    *info = 0;
#line 288 "zlatps.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 289 "zlatps.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 290 "zlatps.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 294 "zlatps.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 295 "zlatps.f"
	*info = -1;
#line 296 "zlatps.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 298 "zlatps.f"
	*info = -2;
#line 299 "zlatps.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 300 "zlatps.f"
	*info = -3;
#line 301 "zlatps.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 303 "zlatps.f"
	*info = -4;
#line 304 "zlatps.f"
    } else if (*n < 0) {
#line 305 "zlatps.f"
	*info = -5;
#line 306 "zlatps.f"
    }
#line 307 "zlatps.f"
    if (*info != 0) {
#line 308 "zlatps.f"
	i__1 = -(*info);
#line 308 "zlatps.f"
	xerbla_("ZLATPS", &i__1, (ftnlen)6);
#line 309 "zlatps.f"
	return 0;
#line 310 "zlatps.f"
    }

/*     Quick return if possible */

#line 314 "zlatps.f"
    if (*n == 0) {
#line 314 "zlatps.f"
	return 0;
#line 314 "zlatps.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 319 "zlatps.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 320 "zlatps.f"
    bignum = 1. / smlnum;
#line 321 "zlatps.f"
    dlabad_(&smlnum, &bignum);
#line 322 "zlatps.f"
    smlnum /= dlamch_("Precision", (ftnlen)9);
#line 323 "zlatps.f"
    bignum = 1. / smlnum;
#line 324 "zlatps.f"
    *scale = 1.;

#line 326 "zlatps.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 330 "zlatps.f"
	if (upper) {

/*           A is upper triangular. */

#line 334 "zlatps.f"
	    ip = 1;
#line 335 "zlatps.f"
	    i__1 = *n;
#line 335 "zlatps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 336 "zlatps.f"
		i__2 = j - 1;
#line 336 "zlatps.f"
		cnorm[j] = dzasum_(&i__2, &ap[ip], &c__1);
#line 337 "zlatps.f"
		ip += j;
#line 338 "zlatps.f"
/* L10: */
#line 338 "zlatps.f"
	    }
#line 339 "zlatps.f"
	} else {

/*           A is lower triangular. */

#line 343 "zlatps.f"
	    ip = 1;
#line 344 "zlatps.f"
	    i__1 = *n - 1;
#line 344 "zlatps.f"
	    for (j = 1; j <= i__1; ++j) {
#line 345 "zlatps.f"
		i__2 = *n - j;
#line 345 "zlatps.f"
		cnorm[j] = dzasum_(&i__2, &ap[ip + 1], &c__1);
#line 346 "zlatps.f"
		ip = ip + *n - j + 1;
#line 347 "zlatps.f"
/* L20: */
#line 347 "zlatps.f"
	    }
#line 348 "zlatps.f"
	    cnorm[*n] = 0.;
#line 349 "zlatps.f"
	}
#line 350 "zlatps.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM/2. */

#line 355 "zlatps.f"
    imax = idamax_(n, &cnorm[1], &c__1);
#line 356 "zlatps.f"
    tmax = cnorm[imax];
#line 357 "zlatps.f"
    if (tmax <= bignum * .5) {
#line 358 "zlatps.f"
	tscal = 1.;
#line 359 "zlatps.f"
    } else {
#line 360 "zlatps.f"
	tscal = .5 / (smlnum * tmax);
#line 361 "zlatps.f"
	dscal_(n, &tscal, &cnorm[1], &c__1);
#line 362 "zlatps.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine ZTPSV can be used. */

#line 367 "zlatps.f"
    xmax = 0.;
#line 368 "zlatps.f"
    i__1 = *n;
#line 368 "zlatps.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 369 "zlatps.f"
	i__2 = j;
#line 369 "zlatps.f"
	d__3 = xmax, d__4 = (d__1 = x[i__2].r / 2., abs(d__1)) + (d__2 = 
		d_imag(&x[j]) / 2., abs(d__2));
#line 369 "zlatps.f"
	xmax = max(d__3,d__4);
#line 370 "zlatps.f"
/* L30: */
#line 370 "zlatps.f"
    }
#line 371 "zlatps.f"
    xbnd = xmax;
#line 372 "zlatps.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 376 "zlatps.f"
	if (upper) {
#line 377 "zlatps.f"
	    jfirst = *n;
#line 378 "zlatps.f"
	    jlast = 1;
#line 379 "zlatps.f"
	    jinc = -1;
#line 380 "zlatps.f"
	} else {
#line 381 "zlatps.f"
	    jfirst = 1;
#line 382 "zlatps.f"
	    jlast = *n;
#line 383 "zlatps.f"
	    jinc = 1;
#line 384 "zlatps.f"
	}

#line 386 "zlatps.f"
	if (tscal != 1.) {
#line 387 "zlatps.f"
	    grow = 0.;
#line 388 "zlatps.f"
	    goto L60;
#line 389 "zlatps.f"
	}

#line 391 "zlatps.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 398 "zlatps.f"
	    grow = .5 / max(xbnd,smlnum);
#line 399 "zlatps.f"
	    xbnd = grow;
#line 400 "zlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 401 "zlatps.f"
	    jlen = *n;
#line 402 "zlatps.f"
	    i__1 = jlast;
#line 402 "zlatps.f"
	    i__2 = jinc;
#line 402 "zlatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 406 "zlatps.f"
		if (grow <= smlnum) {
#line 406 "zlatps.f"
		    goto L60;
#line 406 "zlatps.f"
		}

#line 409 "zlatps.f"
		i__3 = ip;
#line 409 "zlatps.f"
		tjjs.r = ap[i__3].r, tjjs.i = ap[i__3].i;
#line 410 "zlatps.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 412 "zlatps.f"
		if (tjj >= smlnum) {

/*                 M(j) = G(j-1) / abs(A(j,j)) */

/* Computing MIN */
#line 416 "zlatps.f"
		    d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 416 "zlatps.f"
		    xbnd = min(d__1,d__2);
#line 417 "zlatps.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 421 "zlatps.f"
		    xbnd = 0.;
#line 422 "zlatps.f"
		}

#line 424 "zlatps.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 428 "zlatps.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 429 "zlatps.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 433 "zlatps.f"
		    grow = 0.;
#line 434 "zlatps.f"
		}
#line 435 "zlatps.f"
		ip += jinc * jlen;
#line 436 "zlatps.f"
		--jlen;
#line 437 "zlatps.f"
/* L40: */
#line 437 "zlatps.f"
	    }
#line 438 "zlatps.f"
	    grow = xbnd;
#line 439 "zlatps.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 445 "zlatps.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 445 "zlatps.f"
	    grow = min(d__1,d__2);
#line 446 "zlatps.f"
	    i__2 = jlast;
#line 446 "zlatps.f"
	    i__1 = jinc;
#line 446 "zlatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 450 "zlatps.f"
		if (grow <= smlnum) {
#line 450 "zlatps.f"
		    goto L60;
#line 450 "zlatps.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 455 "zlatps.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 456 "zlatps.f"
/* L50: */
#line 456 "zlatps.f"
	    }
#line 457 "zlatps.f"
	}
#line 458 "zlatps.f"
L60:

#line 460 "zlatps.f"
	;
#line 460 "zlatps.f"
    } else {

/*        Compute the growth in A**T * x = b  or  A**H * x = b. */

#line 464 "zlatps.f"
	if (upper) {
#line 465 "zlatps.f"
	    jfirst = 1;
#line 466 "zlatps.f"
	    jlast = *n;
#line 467 "zlatps.f"
	    jinc = 1;
#line 468 "zlatps.f"
	} else {
#line 469 "zlatps.f"
	    jfirst = *n;
#line 470 "zlatps.f"
	    jlast = 1;
#line 471 "zlatps.f"
	    jinc = -1;
#line 472 "zlatps.f"
	}

#line 474 "zlatps.f"
	if (tscal != 1.) {
#line 475 "zlatps.f"
	    grow = 0.;
#line 476 "zlatps.f"
	    goto L90;
#line 477 "zlatps.f"
	}

#line 479 "zlatps.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 486 "zlatps.f"
	    grow = .5 / max(xbnd,smlnum);
#line 487 "zlatps.f"
	    xbnd = grow;
#line 488 "zlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 489 "zlatps.f"
	    jlen = 1;
#line 490 "zlatps.f"
	    i__1 = jlast;
#line 490 "zlatps.f"
	    i__2 = jinc;
#line 490 "zlatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 494 "zlatps.f"
		if (grow <= smlnum) {
#line 494 "zlatps.f"
		    goto L90;
#line 494 "zlatps.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 499 "zlatps.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 500 "zlatps.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 500 "zlatps.f"
		grow = min(d__1,d__2);

#line 502 "zlatps.f"
		i__3 = ip;
#line 502 "zlatps.f"
		tjjs.r = ap[i__3].r, tjjs.i = ap[i__3].i;
#line 503 "zlatps.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 505 "zlatps.f"
		if (tjj >= smlnum) {

/*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 509 "zlatps.f"
		    if (xj > tjj) {
#line 509 "zlatps.f"
			xbnd *= tjj / xj;
#line 509 "zlatps.f"
		    }
#line 511 "zlatps.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 515 "zlatps.f"
		    xbnd = 0.;
#line 516 "zlatps.f"
		}
#line 517 "zlatps.f"
		++jlen;
#line 518 "zlatps.f"
		ip += jinc * jlen;
#line 519 "zlatps.f"
/* L70: */
#line 519 "zlatps.f"
	    }
#line 520 "zlatps.f"
	    grow = min(grow,xbnd);
#line 521 "zlatps.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 527 "zlatps.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 527 "zlatps.f"
	    grow = min(d__1,d__2);
#line 528 "zlatps.f"
	    i__2 = jlast;
#line 528 "zlatps.f"
	    i__1 = jinc;
#line 528 "zlatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 532 "zlatps.f"
		if (grow <= smlnum) {
#line 532 "zlatps.f"
		    goto L90;
#line 532 "zlatps.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 537 "zlatps.f"
		xj = cnorm[j] + 1.;
#line 538 "zlatps.f"
		grow /= xj;
#line 539 "zlatps.f"
/* L80: */
#line 539 "zlatps.f"
	    }
#line 540 "zlatps.f"
	}
#line 541 "zlatps.f"
L90:
#line 542 "zlatps.f"
	;
#line 542 "zlatps.f"
    }

#line 544 "zlatps.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 549 "zlatps.f"
	ztpsv_(uplo, trans, diag, n, &ap[1], &x[1], &c__1, (ftnlen)1, (ftnlen)
		1, (ftnlen)1);
#line 550 "zlatps.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 554 "zlatps.f"
	if (xmax > bignum * .5) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 559 "zlatps.f"
	    *scale = bignum * .5 / xmax;
#line 560 "zlatps.f"
	    zdscal_(n, scale, &x[1], &c__1);
#line 561 "zlatps.f"
	    xmax = bignum;
#line 562 "zlatps.f"
	} else {
#line 563 "zlatps.f"
	    xmax *= 2.;
#line 564 "zlatps.f"
	}

#line 566 "zlatps.f"
	if (notran) {

/*           Solve A * x = b */

#line 570 "zlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 571 "zlatps.f"
	    i__1 = jlast;
#line 571 "zlatps.f"
	    i__2 = jinc;
#line 571 "zlatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 575 "zlatps.f"
		i__3 = j;
#line 575 "zlatps.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 576 "zlatps.f"
		if (nounit) {
#line 577 "zlatps.f"
		    i__3 = ip;
#line 577 "zlatps.f"
		    z__1.r = tscal * ap[i__3].r, z__1.i = tscal * ap[i__3].i;
#line 577 "zlatps.f"
		    tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 578 "zlatps.f"
		} else {
#line 579 "zlatps.f"
		    tjjs.r = tscal, tjjs.i = 0.;
#line 580 "zlatps.f"
		    if (tscal == 1.) {
#line 580 "zlatps.f"
			goto L110;
#line 580 "zlatps.f"
		    }
#line 582 "zlatps.f"
		}
#line 583 "zlatps.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));
#line 584 "zlatps.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 588 "zlatps.f"
		    if (tjj < 1.) {
#line 589 "zlatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 593 "zlatps.f"
			    rec = 1. / xj;
#line 594 "zlatps.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 595 "zlatps.f"
			    *scale *= rec;
#line 596 "zlatps.f"
			    xmax *= rec;
#line 597 "zlatps.f"
			}
#line 598 "zlatps.f"
		    }
#line 599 "zlatps.f"
		    i__3 = j;
#line 599 "zlatps.f"
		    zladiv_(&z__1, &x[j], &tjjs);
#line 599 "zlatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 600 "zlatps.f"
		    i__3 = j;
#line 600 "zlatps.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 601 "zlatps.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 605 "zlatps.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 610 "zlatps.f"
			rec = tjj * bignum / xj;
#line 611 "zlatps.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 616 "zlatps.f"
			    rec /= cnorm[j];
#line 617 "zlatps.f"
			}
#line 618 "zlatps.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 619 "zlatps.f"
			*scale *= rec;
#line 620 "zlatps.f"
			xmax *= rec;
#line 621 "zlatps.f"
		    }
#line 622 "zlatps.f"
		    i__3 = j;
#line 622 "zlatps.f"
		    zladiv_(&z__1, &x[j], &tjjs);
#line 622 "zlatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 623 "zlatps.f"
		    i__3 = j;
#line 623 "zlatps.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 624 "zlatps.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 629 "zlatps.f"
		    i__3 = *n;
#line 629 "zlatps.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 630 "zlatps.f"
			i__4 = i__;
#line 630 "zlatps.f"
			x[i__4].r = 0., x[i__4].i = 0.;
#line 631 "zlatps.f"
/* L100: */
#line 631 "zlatps.f"
		    }
#line 632 "zlatps.f"
		    i__3 = j;
#line 632 "zlatps.f"
		    x[i__3].r = 1., x[i__3].i = 0.;
#line 633 "zlatps.f"
		    xj = 1.;
#line 634 "zlatps.f"
		    *scale = 0.;
#line 635 "zlatps.f"
		    xmax = 0.;
#line 636 "zlatps.f"
		}
#line 637 "zlatps.f"
L110:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 642 "zlatps.f"
		if (xj > 1.) {
#line 643 "zlatps.f"
		    rec = 1. / xj;
#line 644 "zlatps.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 648 "zlatps.f"
			rec *= .5;
#line 649 "zlatps.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 650 "zlatps.f"
			*scale *= rec;
#line 651 "zlatps.f"
		    }
#line 652 "zlatps.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 656 "zlatps.f"
		    zdscal_(n, &c_b36, &x[1], &c__1);
#line 657 "zlatps.f"
		    *scale *= .5;
#line 658 "zlatps.f"
		}

#line 660 "zlatps.f"
		if (upper) {
#line 661 "zlatps.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

#line 666 "zlatps.f"
			i__3 = j - 1;
#line 666 "zlatps.f"
			i__4 = j;
#line 666 "zlatps.f"
			z__2.r = -x[i__4].r, z__2.i = -x[i__4].i;
#line 666 "zlatps.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 666 "zlatps.f"
			zaxpy_(&i__3, &z__1, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 668 "zlatps.f"
			i__3 = j - 1;
#line 668 "zlatps.f"
			i__ = izamax_(&i__3, &x[1], &c__1);
#line 669 "zlatps.f"
			i__3 = i__;
#line 669 "zlatps.f"
			xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&x[i__]), abs(d__2));
#line 670 "zlatps.f"
		    }
#line 671 "zlatps.f"
		    ip -= j;
#line 672 "zlatps.f"
		} else {
#line 673 "zlatps.f"
		    if (j < *n) {

/*                    Compute the update */
/*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

#line 678 "zlatps.f"
			i__3 = *n - j;
#line 678 "zlatps.f"
			i__4 = j;
#line 678 "zlatps.f"
			z__2.r = -x[i__4].r, z__2.i = -x[i__4].i;
#line 678 "zlatps.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 678 "zlatps.f"
			zaxpy_(&i__3, &z__1, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 680 "zlatps.f"
			i__3 = *n - j;
#line 680 "zlatps.f"
			i__ = j + izamax_(&i__3, &x[j + 1], &c__1);
#line 681 "zlatps.f"
			i__3 = i__;
#line 681 "zlatps.f"
			xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&x[i__]), abs(d__2));
#line 682 "zlatps.f"
		    }
#line 683 "zlatps.f"
		    ip = ip + *n - j + 1;
#line 684 "zlatps.f"
		}
#line 685 "zlatps.f"
/* L120: */
#line 685 "zlatps.f"
	    }

#line 687 "zlatps.f"
	} else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*           Solve A**T * x = b */

#line 691 "zlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 692 "zlatps.f"
	    jlen = 1;
#line 693 "zlatps.f"
	    i__2 = jlast;
#line 693 "zlatps.f"
	    i__1 = jinc;
#line 693 "zlatps.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 698 "zlatps.f"
		i__3 = j;
#line 698 "zlatps.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 699 "zlatps.f"
		uscal.r = tscal, uscal.i = 0.;
#line 700 "zlatps.f"
		rec = 1. / max(xmax,1.);
#line 701 "zlatps.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 705 "zlatps.f"
		    rec *= .5;
#line 706 "zlatps.f"
		    if (nounit) {
#line 707 "zlatps.f"
			i__3 = ip;
#line 707 "zlatps.f"
			z__1.r = tscal * ap[i__3].r, z__1.i = tscal * ap[i__3]
				.i;
#line 707 "zlatps.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 708 "zlatps.f"
		    } else {
#line 709 "zlatps.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 710 "zlatps.f"
		    }
#line 711 "zlatps.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 712 "zlatps.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 716 "zlatps.f"
			d__1 = 1., d__2 = rec * tjj;
#line 716 "zlatps.f"
			rec = min(d__1,d__2);
#line 717 "zlatps.f"
			zladiv_(&z__1, &uscal, &tjjs);
#line 717 "zlatps.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 718 "zlatps.f"
		    }
#line 719 "zlatps.f"
		    if (rec < 1.) {
#line 720 "zlatps.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 721 "zlatps.f"
			*scale *= rec;
#line 722 "zlatps.f"
			xmax *= rec;
#line 723 "zlatps.f"
		    }
#line 724 "zlatps.f"
		}

#line 726 "zlatps.f"
		csumj.r = 0., csumj.i = 0.;
#line 727 "zlatps.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call ZDOTU to perform the dot product. */

#line 732 "zlatps.f"
		    if (upper) {
#line 733 "zlatps.f"
			i__3 = j - 1;
#line 733 "zlatps.f"
			zdotu_(&z__1, &i__3, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 733 "zlatps.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 734 "zlatps.f"
		    } else if (j < *n) {
#line 735 "zlatps.f"
			i__3 = *n - j;
#line 735 "zlatps.f"
			zdotu_(&z__1, &i__3, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 735 "zlatps.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 736 "zlatps.f"
		    }
#line 737 "zlatps.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 741 "zlatps.f"
		    if (upper) {
#line 742 "zlatps.f"
			i__3 = j - 1;
#line 742 "zlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 743 "zlatps.f"
			    i__4 = ip - j + i__;
#line 743 "zlatps.f"
			    z__3.r = ap[i__4].r * uscal.r - ap[i__4].i * 
				    uscal.i, z__3.i = ap[i__4].r * uscal.i + 
				    ap[i__4].i * uscal.r;
#line 743 "zlatps.f"
			    i__5 = i__;
#line 743 "zlatps.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 743 "zlatps.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 743 "zlatps.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 744 "zlatps.f"
/* L130: */
#line 744 "zlatps.f"
			}
#line 745 "zlatps.f"
		    } else if (j < *n) {
#line 746 "zlatps.f"
			i__3 = *n - j;
#line 746 "zlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 747 "zlatps.f"
			    i__4 = ip + i__;
#line 747 "zlatps.f"
			    z__3.r = ap[i__4].r * uscal.r - ap[i__4].i * 
				    uscal.i, z__3.i = ap[i__4].r * uscal.i + 
				    ap[i__4].i * uscal.r;
#line 747 "zlatps.f"
			    i__5 = j + i__;
#line 747 "zlatps.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 747 "zlatps.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 747 "zlatps.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 748 "zlatps.f"
/* L140: */
#line 748 "zlatps.f"
			}
#line 749 "zlatps.f"
		    }
#line 750 "zlatps.f"
		}

#line 752 "zlatps.f"
		z__1.r = tscal, z__1.i = 0.;
#line 752 "zlatps.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 757 "zlatps.f"
		    i__3 = j;
#line 757 "zlatps.f"
		    i__4 = j;
#line 757 "zlatps.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 757 "zlatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 758 "zlatps.f"
		    i__3 = j;
#line 758 "zlatps.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 759 "zlatps.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 763 "zlatps.f"
			i__3 = ip;
#line 763 "zlatps.f"
			z__1.r = tscal * ap[i__3].r, z__1.i = tscal * ap[i__3]
				.i;
#line 763 "zlatps.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 764 "zlatps.f"
		    } else {
#line 765 "zlatps.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 766 "zlatps.f"
			if (tscal == 1.) {
#line 766 "zlatps.f"
			    goto L160;
#line 766 "zlatps.f"
			}
#line 768 "zlatps.f"
		    }
#line 769 "zlatps.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 770 "zlatps.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 774 "zlatps.f"
			if (tjj < 1.) {
#line 775 "zlatps.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 779 "zlatps.f"
				rec = 1. / xj;
#line 780 "zlatps.f"
				zdscal_(n, &rec, &x[1], &c__1);
#line 781 "zlatps.f"
				*scale *= rec;
#line 782 "zlatps.f"
				xmax *= rec;
#line 783 "zlatps.f"
			    }
#line 784 "zlatps.f"
			}
#line 785 "zlatps.f"
			i__3 = j;
#line 785 "zlatps.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 785 "zlatps.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 786 "zlatps.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 790 "zlatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 794 "zlatps.f"
			    rec = tjj * bignum / xj;
#line 795 "zlatps.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 796 "zlatps.f"
			    *scale *= rec;
#line 797 "zlatps.f"
			    xmax *= rec;
#line 798 "zlatps.f"
			}
#line 799 "zlatps.f"
			i__3 = j;
#line 799 "zlatps.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 799 "zlatps.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 800 "zlatps.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**T *x = 0. */

#line 805 "zlatps.f"
			i__3 = *n;
#line 805 "zlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 806 "zlatps.f"
			    i__4 = i__;
#line 806 "zlatps.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 807 "zlatps.f"
/* L150: */
#line 807 "zlatps.f"
			}
#line 808 "zlatps.f"
			i__3 = j;
#line 808 "zlatps.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 809 "zlatps.f"
			*scale = 0.;
#line 810 "zlatps.f"
			xmax = 0.;
#line 811 "zlatps.f"
		    }
#line 812 "zlatps.f"
L160:
#line 813 "zlatps.f"
		    ;
#line 813 "zlatps.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 818 "zlatps.f"
		    i__3 = j;
#line 818 "zlatps.f"
		    zladiv_(&z__2, &x[j], &tjjs);
#line 818 "zlatps.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 818 "zlatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 819 "zlatps.f"
		}
/* Computing MAX */
#line 820 "zlatps.f"
		i__3 = j;
#line 820 "zlatps.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 820 "zlatps.f"
		xmax = max(d__3,d__4);
#line 821 "zlatps.f"
		++jlen;
#line 822 "zlatps.f"
		ip += jinc * jlen;
#line 823 "zlatps.f"
/* L170: */
#line 823 "zlatps.f"
	    }

#line 825 "zlatps.f"
	} else {

/*           Solve A**H * x = b */

#line 829 "zlatps.f"
	    ip = jfirst * (jfirst + 1) / 2;
#line 830 "zlatps.f"
	    jlen = 1;
#line 831 "zlatps.f"
	    i__1 = jlast;
#line 831 "zlatps.f"
	    i__2 = jinc;
#line 831 "zlatps.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 836 "zlatps.f"
		i__3 = j;
#line 836 "zlatps.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 837 "zlatps.f"
		uscal.r = tscal, uscal.i = 0.;
#line 838 "zlatps.f"
		rec = 1. / max(xmax,1.);
#line 839 "zlatps.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 843 "zlatps.f"
		    rec *= .5;
#line 844 "zlatps.f"
		    if (nounit) {
#line 845 "zlatps.f"
			d_cnjg(&z__2, &ap[ip]);
#line 845 "zlatps.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 845 "zlatps.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 846 "zlatps.f"
		    } else {
#line 847 "zlatps.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 848 "zlatps.f"
		    }
#line 849 "zlatps.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 850 "zlatps.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 854 "zlatps.f"
			d__1 = 1., d__2 = rec * tjj;
#line 854 "zlatps.f"
			rec = min(d__1,d__2);
#line 855 "zlatps.f"
			zladiv_(&z__1, &uscal, &tjjs);
#line 855 "zlatps.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 856 "zlatps.f"
		    }
#line 857 "zlatps.f"
		    if (rec < 1.) {
#line 858 "zlatps.f"
			zdscal_(n, &rec, &x[1], &c__1);
#line 859 "zlatps.f"
			*scale *= rec;
#line 860 "zlatps.f"
			xmax *= rec;
#line 861 "zlatps.f"
		    }
#line 862 "zlatps.f"
		}

#line 864 "zlatps.f"
		csumj.r = 0., csumj.i = 0.;
#line 865 "zlatps.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call ZDOTC to perform the dot product. */

#line 870 "zlatps.f"
		    if (upper) {
#line 871 "zlatps.f"
			i__3 = j - 1;
#line 871 "zlatps.f"
			zdotc_(&z__1, &i__3, &ap[ip - j + 1], &c__1, &x[1], &
				c__1);
#line 871 "zlatps.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 872 "zlatps.f"
		    } else if (j < *n) {
#line 873 "zlatps.f"
			i__3 = *n - j;
#line 873 "zlatps.f"
			zdotc_(&z__1, &i__3, &ap[ip + 1], &c__1, &x[j + 1], &
				c__1);
#line 873 "zlatps.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 874 "zlatps.f"
		    }
#line 875 "zlatps.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 879 "zlatps.f"
		    if (upper) {
#line 880 "zlatps.f"
			i__3 = j - 1;
#line 880 "zlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 881 "zlatps.f"
			    d_cnjg(&z__4, &ap[ip - j + i__]);
#line 881 "zlatps.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 881 "zlatps.f"
			    i__4 = i__;
#line 881 "zlatps.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 881 "zlatps.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 881 "zlatps.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 883 "zlatps.f"
/* L180: */
#line 883 "zlatps.f"
			}
#line 884 "zlatps.f"
		    } else if (j < *n) {
#line 885 "zlatps.f"
			i__3 = *n - j;
#line 885 "zlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 886 "zlatps.f"
			    d_cnjg(&z__4, &ap[ip + i__]);
#line 886 "zlatps.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 886 "zlatps.f"
			    i__4 = j + i__;
#line 886 "zlatps.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 886 "zlatps.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 886 "zlatps.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 888 "zlatps.f"
/* L190: */
#line 888 "zlatps.f"
			}
#line 889 "zlatps.f"
		    }
#line 890 "zlatps.f"
		}

#line 892 "zlatps.f"
		z__1.r = tscal, z__1.i = 0.;
#line 892 "zlatps.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 897 "zlatps.f"
		    i__3 = j;
#line 897 "zlatps.f"
		    i__4 = j;
#line 897 "zlatps.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 897 "zlatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 898 "zlatps.f"
		    i__3 = j;
#line 898 "zlatps.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 899 "zlatps.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 903 "zlatps.f"
			d_cnjg(&z__2, &ap[ip]);
#line 903 "zlatps.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 903 "zlatps.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 904 "zlatps.f"
		    } else {
#line 905 "zlatps.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 906 "zlatps.f"
			if (tscal == 1.) {
#line 906 "zlatps.f"
			    goto L210;
#line 906 "zlatps.f"
			}
#line 908 "zlatps.f"
		    }
#line 909 "zlatps.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 910 "zlatps.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 914 "zlatps.f"
			if (tjj < 1.) {
#line 915 "zlatps.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 919 "zlatps.f"
				rec = 1. / xj;
#line 920 "zlatps.f"
				zdscal_(n, &rec, &x[1], &c__1);
#line 921 "zlatps.f"
				*scale *= rec;
#line 922 "zlatps.f"
				xmax *= rec;
#line 923 "zlatps.f"
			    }
#line 924 "zlatps.f"
			}
#line 925 "zlatps.f"
			i__3 = j;
#line 925 "zlatps.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 925 "zlatps.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 926 "zlatps.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 930 "zlatps.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 934 "zlatps.f"
			    rec = tjj * bignum / xj;
#line 935 "zlatps.f"
			    zdscal_(n, &rec, &x[1], &c__1);
#line 936 "zlatps.f"
			    *scale *= rec;
#line 937 "zlatps.f"
			    xmax *= rec;
#line 938 "zlatps.f"
			}
#line 939 "zlatps.f"
			i__3 = j;
#line 939 "zlatps.f"
			zladiv_(&z__1, &x[j], &tjjs);
#line 939 "zlatps.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 940 "zlatps.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**H *x = 0. */

#line 945 "zlatps.f"
			i__3 = *n;
#line 945 "zlatps.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 946 "zlatps.f"
			    i__4 = i__;
#line 946 "zlatps.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 947 "zlatps.f"
/* L200: */
#line 947 "zlatps.f"
			}
#line 948 "zlatps.f"
			i__3 = j;
#line 948 "zlatps.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 949 "zlatps.f"
			*scale = 0.;
#line 950 "zlatps.f"
			xmax = 0.;
#line 951 "zlatps.f"
		    }
#line 952 "zlatps.f"
L210:
#line 953 "zlatps.f"
		    ;
#line 953 "zlatps.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 958 "zlatps.f"
		    i__3 = j;
#line 958 "zlatps.f"
		    zladiv_(&z__2, &x[j], &tjjs);
#line 958 "zlatps.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 958 "zlatps.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 959 "zlatps.f"
		}
/* Computing MAX */
#line 960 "zlatps.f"
		i__3 = j;
#line 960 "zlatps.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 960 "zlatps.f"
		xmax = max(d__3,d__4);
#line 961 "zlatps.f"
		++jlen;
#line 962 "zlatps.f"
		ip += jinc * jlen;
#line 963 "zlatps.f"
/* L220: */
#line 963 "zlatps.f"
	    }
#line 964 "zlatps.f"
	}
#line 965 "zlatps.f"
	*scale /= tscal;
#line 966 "zlatps.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 970 "zlatps.f"
    if (tscal != 1.) {
#line 971 "zlatps.f"
	d__1 = 1. / tscal;
#line 971 "zlatps.f"
	dscal_(n, &d__1, &cnorm[1], &c__1);
#line 972 "zlatps.f"
    }

#line 974 "zlatps.f"
    return 0;

/*     End of ZLATPS */

} /* zlatps_ */

