#line 1 "clatbs.f"
/* clatbs.f -- translated by f2c (version 20100827).
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

#line 1 "clatbs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b36 = .5;

/* > \brief \b CLATBS solves a triangular banded system of equations. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLATBS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clatbs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clatbs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clatbs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X, */
/*                          SCALE, CNORM, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, NORMIN, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       REAL               SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               CNORM( * ) */
/*       COMPLEX            AB( LDAB, * ), X( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLATBS solves one of the triangular systems */
/* > */
/* >    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b, */
/* > */
/* > with scaling to prevent overflow, where A is an upper or lower */
/* > triangular band matrix.  Here A**T denotes the transpose of A, x and b */
/* > are n-element vectors, and s is a scaling factor, usually less than */
/* > or equal to 1, chosen so that the components of x will be less than */
/* > the overflow threshold.  If the unscaled problem will not cause */
/* > overflow, the Level 2 BLAS routine CTBSV is called.  If the matrix A */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of subdiagonals or superdiagonals in the */
/* >          triangular matrix A.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          X is COMPLEX array, dimension (N) */
/* >          On entry, the right hand side b of the triangular system. */
/* >          On exit, X is overwritten by the solution vector x. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >          The scaling factor s for the triangular system */
/* >             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b. */
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

/* > \ingroup complexOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A rough bound on x is computed; if that is less than overflow, CTBSV */
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
/* >  Since |x(j)| <= M(j), we use the Level 2 BLAS routine CTBSV if the */
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
/* >  and we can safely call CTBSV if 1/M(n) and 1/G(n) are both greater */
/* >  than max(underflow, 1/overflow). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int clatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublecomplex *ab, integer *ldab, 
	doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info,
	 ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len, ftnlen 
	normin_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal xj, rec, tjj;
    static integer jinc, jlen;
    static doublereal xbnd;
    static integer imax;
    static doublereal tmax;
    static doublecomplex tjjs;
    static doublereal xmax, grow;
    static integer maind;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal tscal;
    static doublecomplex uscal;
    static integer jlast;
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublecomplex csumj;
    extern /* Subroutine */ int ctbsv_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern /* Double Complex */ VOID cladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 299 "clatbs.f"
    /* Parameter adjustments */
#line 299 "clatbs.f"
    ab_dim1 = *ldab;
#line 299 "clatbs.f"
    ab_offset = 1 + ab_dim1;
#line 299 "clatbs.f"
    ab -= ab_offset;
#line 299 "clatbs.f"
    --x;
#line 299 "clatbs.f"
    --cnorm;
#line 299 "clatbs.f"

#line 299 "clatbs.f"
    /* Function Body */
#line 299 "clatbs.f"
    *info = 0;
#line 300 "clatbs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 301 "clatbs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 302 "clatbs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 306 "clatbs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 307 "clatbs.f"
	*info = -1;
#line 308 "clatbs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 310 "clatbs.f"
	*info = -2;
#line 311 "clatbs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 312 "clatbs.f"
	*info = -3;
#line 313 "clatbs.f"
    } else if (! lsame_(normin, "Y", (ftnlen)1, (ftnlen)1) && ! lsame_(normin,
	     "N", (ftnlen)1, (ftnlen)1)) {
#line 315 "clatbs.f"
	*info = -4;
#line 316 "clatbs.f"
    } else if (*n < 0) {
#line 317 "clatbs.f"
	*info = -5;
#line 318 "clatbs.f"
    } else if (*kd < 0) {
#line 319 "clatbs.f"
	*info = -6;
#line 320 "clatbs.f"
    } else if (*ldab < *kd + 1) {
#line 321 "clatbs.f"
	*info = -8;
#line 322 "clatbs.f"
    }
#line 323 "clatbs.f"
    if (*info != 0) {
#line 324 "clatbs.f"
	i__1 = -(*info);
#line 324 "clatbs.f"
	xerbla_("CLATBS", &i__1, (ftnlen)6);
#line 325 "clatbs.f"
	return 0;
#line 326 "clatbs.f"
    }

/*     Quick return if possible */

#line 330 "clatbs.f"
    if (*n == 0) {
#line 330 "clatbs.f"
	return 0;
#line 330 "clatbs.f"
    }

/*     Determine machine dependent parameters to control overflow. */

#line 335 "clatbs.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 336 "clatbs.f"
    bignum = 1. / smlnum;
#line 337 "clatbs.f"
    slabad_(&smlnum, &bignum);
#line 338 "clatbs.f"
    smlnum /= slamch_("Precision", (ftnlen)9);
#line 339 "clatbs.f"
    bignum = 1. / smlnum;
#line 340 "clatbs.f"
    *scale = 1.;

#line 342 "clatbs.f"
    if (lsame_(normin, "N", (ftnlen)1, (ftnlen)1)) {

/*        Compute the 1-norm of each column, not including the diagonal. */

#line 346 "clatbs.f"
	if (upper) {

/*           A is upper triangular. */

#line 350 "clatbs.f"
	    i__1 = *n;
#line 350 "clatbs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 351 "clatbs.f"
		i__2 = *kd, i__3 = j - 1;
#line 351 "clatbs.f"
		jlen = min(i__2,i__3);
#line 352 "clatbs.f"
		cnorm[j] = scasum_(&jlen, &ab[*kd + 1 - jlen + j * ab_dim1], &
			c__1);
#line 353 "clatbs.f"
/* L10: */
#line 353 "clatbs.f"
	    }
#line 354 "clatbs.f"
	} else {

/*           A is lower triangular. */

#line 358 "clatbs.f"
	    i__1 = *n;
#line 358 "clatbs.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 359 "clatbs.f"
		i__2 = *kd, i__3 = *n - j;
#line 359 "clatbs.f"
		jlen = min(i__2,i__3);
#line 360 "clatbs.f"
		if (jlen > 0) {
#line 361 "clatbs.f"
		    cnorm[j] = scasum_(&jlen, &ab[j * ab_dim1 + 2], &c__1);
#line 362 "clatbs.f"
		} else {
#line 363 "clatbs.f"
		    cnorm[j] = 0.;
#line 364 "clatbs.f"
		}
#line 365 "clatbs.f"
/* L20: */
#line 365 "clatbs.f"
	    }
#line 366 "clatbs.f"
	}
#line 367 "clatbs.f"
    }

/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
/*     greater than BIGNUM/2. */

#line 372 "clatbs.f"
    imax = isamax_(n, &cnorm[1], &c__1);
#line 373 "clatbs.f"
    tmax = cnorm[imax];
#line 374 "clatbs.f"
    if (tmax <= bignum * .5) {
#line 375 "clatbs.f"
	tscal = 1.;
#line 376 "clatbs.f"
    } else {
#line 377 "clatbs.f"
	tscal = .5 / (smlnum * tmax);
#line 378 "clatbs.f"
	sscal_(n, &tscal, &cnorm[1], &c__1);
#line 379 "clatbs.f"
    }

/*     Compute a bound on the computed solution vector to see if the */
/*     Level 2 BLAS routine CTBSV can be used. */

#line 384 "clatbs.f"
    xmax = 0.;
#line 385 "clatbs.f"
    i__1 = *n;
#line 385 "clatbs.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 386 "clatbs.f"
	i__2 = j;
#line 386 "clatbs.f"
	d__3 = xmax, d__4 = (d__1 = x[i__2].r / 2., abs(d__1)) + (d__2 = 
		d_imag(&x[j]) / 2., abs(d__2));
#line 386 "clatbs.f"
	xmax = max(d__3,d__4);
#line 387 "clatbs.f"
/* L30: */
#line 387 "clatbs.f"
    }
#line 388 "clatbs.f"
    xbnd = xmax;
#line 389 "clatbs.f"
    if (notran) {

/*        Compute the growth in A * x = b. */

#line 393 "clatbs.f"
	if (upper) {
#line 394 "clatbs.f"
	    jfirst = *n;
#line 395 "clatbs.f"
	    jlast = 1;
#line 396 "clatbs.f"
	    jinc = -1;
#line 397 "clatbs.f"
	    maind = *kd + 1;
#line 398 "clatbs.f"
	} else {
#line 399 "clatbs.f"
	    jfirst = 1;
#line 400 "clatbs.f"
	    jlast = *n;
#line 401 "clatbs.f"
	    jinc = 1;
#line 402 "clatbs.f"
	    maind = 1;
#line 403 "clatbs.f"
	}

#line 405 "clatbs.f"
	if (tscal != 1.) {
#line 406 "clatbs.f"
	    grow = 0.;
#line 407 "clatbs.f"
	    goto L60;
#line 408 "clatbs.f"
	}

#line 410 "clatbs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, G(0) = max{x(i), i=1,...,n}. */

#line 417 "clatbs.f"
	    grow = .5 / max(xbnd,smlnum);
#line 418 "clatbs.f"
	    xbnd = grow;
#line 419 "clatbs.f"
	    i__1 = jlast;
#line 419 "clatbs.f"
	    i__2 = jinc;
#line 419 "clatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 423 "clatbs.f"
		if (grow <= smlnum) {
#line 423 "clatbs.f"
		    goto L60;
#line 423 "clatbs.f"
		}

#line 426 "clatbs.f"
		i__3 = maind + j * ab_dim1;
#line 426 "clatbs.f"
		tjjs.r = ab[i__3].r, tjjs.i = ab[i__3].i;
#line 427 "clatbs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 429 "clatbs.f"
		if (tjj >= smlnum) {

/*                 M(j) = G(j-1) / abs(A(j,j)) */

/* Computing MIN */
#line 433 "clatbs.f"
		    d__1 = xbnd, d__2 = min(1.,tjj) * grow;
#line 433 "clatbs.f"
		    xbnd = min(d__1,d__2);
#line 434 "clatbs.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 438 "clatbs.f"
		    xbnd = 0.;
#line 439 "clatbs.f"
		}

#line 441 "clatbs.f"
		if (tjj + cnorm[j] >= smlnum) {

/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

#line 445 "clatbs.f"
		    grow *= tjj / (tjj + cnorm[j]);
#line 446 "clatbs.f"
		} else {

/*                 G(j) could overflow, set GROW to 0. */

#line 450 "clatbs.f"
		    grow = 0.;
#line 451 "clatbs.f"
		}
#line 452 "clatbs.f"
/* L40: */
#line 452 "clatbs.f"
	    }
#line 453 "clatbs.f"
	    grow = xbnd;
#line 454 "clatbs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 460 "clatbs.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 460 "clatbs.f"
	    grow = min(d__1,d__2);
#line 461 "clatbs.f"
	    i__2 = jlast;
#line 461 "clatbs.f"
	    i__1 = jinc;
#line 461 "clatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 465 "clatbs.f"
		if (grow <= smlnum) {
#line 465 "clatbs.f"
		    goto L60;
#line 465 "clatbs.f"
		}

/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

#line 470 "clatbs.f"
		grow *= 1. / (cnorm[j] + 1.);
#line 471 "clatbs.f"
/* L50: */
#line 471 "clatbs.f"
	    }
#line 472 "clatbs.f"
	}
#line 473 "clatbs.f"
L60:

#line 475 "clatbs.f"
	;
#line 475 "clatbs.f"
    } else {

/*        Compute the growth in A**T * x = b  or  A**H * x = b. */

#line 479 "clatbs.f"
	if (upper) {
#line 480 "clatbs.f"
	    jfirst = 1;
#line 481 "clatbs.f"
	    jlast = *n;
#line 482 "clatbs.f"
	    jinc = 1;
#line 483 "clatbs.f"
	    maind = *kd + 1;
#line 484 "clatbs.f"
	} else {
#line 485 "clatbs.f"
	    jfirst = *n;
#line 486 "clatbs.f"
	    jlast = 1;
#line 487 "clatbs.f"
	    jinc = -1;
#line 488 "clatbs.f"
	    maind = 1;
#line 489 "clatbs.f"
	}

#line 491 "clatbs.f"
	if (tscal != 1.) {
#line 492 "clatbs.f"
	    grow = 0.;
#line 493 "clatbs.f"
	    goto L90;
#line 494 "clatbs.f"
	}

#line 496 "clatbs.f"
	if (nounit) {

/*           A is non-unit triangular. */

/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
/*           Initially, M(0) = max{x(i), i=1,...,n}. */

#line 503 "clatbs.f"
	    grow = .5 / max(xbnd,smlnum);
#line 504 "clatbs.f"
	    xbnd = grow;
#line 505 "clatbs.f"
	    i__1 = jlast;
#line 505 "clatbs.f"
	    i__2 = jinc;
#line 505 "clatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Exit the loop if the growth factor is too small. */

#line 509 "clatbs.f"
		if (grow <= smlnum) {
#line 509 "clatbs.f"
		    goto L90;
#line 509 "clatbs.f"
		}

/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

#line 514 "clatbs.f"
		xj = cnorm[j] + 1.;
/* Computing MIN */
#line 515 "clatbs.f"
		d__1 = grow, d__2 = xbnd / xj;
#line 515 "clatbs.f"
		grow = min(d__1,d__2);

#line 517 "clatbs.f"
		i__3 = maind + j * ab_dim1;
#line 517 "clatbs.f"
		tjjs.r = ab[i__3].r, tjjs.i = ab[i__3].i;
#line 518 "clatbs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));

#line 520 "clatbs.f"
		if (tjj >= smlnum) {

/*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

#line 524 "clatbs.f"
		    if (xj > tjj) {
#line 524 "clatbs.f"
			xbnd *= tjj / xj;
#line 524 "clatbs.f"
		    }
#line 526 "clatbs.f"
		} else {

/*                 M(j) could overflow, set XBND to 0. */

#line 530 "clatbs.f"
		    xbnd = 0.;
#line 531 "clatbs.f"
		}
#line 532 "clatbs.f"
/* L70: */
#line 532 "clatbs.f"
	    }
#line 533 "clatbs.f"
	    grow = min(grow,xbnd);
#line 534 "clatbs.f"
	} else {

/*           A is unit triangular. */

/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

/* Computing MIN */
#line 540 "clatbs.f"
	    d__1 = 1., d__2 = .5 / max(xbnd,smlnum);
#line 540 "clatbs.f"
	    grow = min(d__1,d__2);
#line 541 "clatbs.f"
	    i__2 = jlast;
#line 541 "clatbs.f"
	    i__1 = jinc;
#line 541 "clatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Exit the loop if the growth factor is too small. */

#line 545 "clatbs.f"
		if (grow <= smlnum) {
#line 545 "clatbs.f"
		    goto L90;
#line 545 "clatbs.f"
		}

/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

#line 550 "clatbs.f"
		xj = cnorm[j] + 1.;
#line 551 "clatbs.f"
		grow /= xj;
#line 552 "clatbs.f"
/* L80: */
#line 552 "clatbs.f"
	    }
#line 553 "clatbs.f"
	}
#line 554 "clatbs.f"
L90:
#line 555 "clatbs.f"
	;
#line 555 "clatbs.f"
    }

#line 557 "clatbs.f"
    if (grow * tscal > smlnum) {

/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
/*        elements of X is not too small. */

#line 562 "clatbs.f"
	ctbsv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &x[1], &c__1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 563 "clatbs.f"
    } else {

/*        Use a Level 1 BLAS solve, scaling intermediate results. */

#line 567 "clatbs.f"
	if (xmax > bignum * .5) {

/*           Scale X so that its components are less than or equal to */
/*           BIGNUM in absolute value. */

#line 572 "clatbs.f"
	    *scale = bignum * .5 / xmax;
#line 573 "clatbs.f"
	    csscal_(n, scale, &x[1], &c__1);
#line 574 "clatbs.f"
	    xmax = bignum;
#line 575 "clatbs.f"
	} else {
#line 576 "clatbs.f"
	    xmax *= 2.;
#line 577 "clatbs.f"
	}

#line 579 "clatbs.f"
	if (notran) {

/*           Solve A * x = b */

#line 583 "clatbs.f"
	    i__1 = jlast;
#line 583 "clatbs.f"
	    i__2 = jinc;
#line 583 "clatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

#line 587 "clatbs.f"
		i__3 = j;
#line 587 "clatbs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 588 "clatbs.f"
		if (nounit) {
#line 589 "clatbs.f"
		    i__3 = maind + j * ab_dim1;
#line 589 "clatbs.f"
		    z__1.r = tscal * ab[i__3].r, z__1.i = tscal * ab[i__3].i;
#line 589 "clatbs.f"
		    tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 590 "clatbs.f"
		} else {
#line 591 "clatbs.f"
		    tjjs.r = tscal, tjjs.i = 0.;
#line 592 "clatbs.f"
		    if (tscal == 1.) {
#line 592 "clatbs.f"
			goto L105;
#line 592 "clatbs.f"
		    }
#line 594 "clatbs.f"
		}
#line 595 "clatbs.f"
		tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), abs(
			d__2));
#line 596 "clatbs.f"
		if (tjj > smlnum) {

/*                    abs(A(j,j)) > SMLNUM: */

#line 600 "clatbs.f"
		    if (tjj < 1.) {
#line 601 "clatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by 1/b(j). */

#line 605 "clatbs.f"
			    rec = 1. / xj;
#line 606 "clatbs.f"
			    csscal_(n, &rec, &x[1], &c__1);
#line 607 "clatbs.f"
			    *scale *= rec;
#line 608 "clatbs.f"
			    xmax *= rec;
#line 609 "clatbs.f"
			}
#line 610 "clatbs.f"
		    }
#line 611 "clatbs.f"
		    i__3 = j;
#line 611 "clatbs.f"
		    cladiv_(&z__1, &x[j], &tjjs);
#line 611 "clatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 612 "clatbs.f"
		    i__3 = j;
#line 612 "clatbs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 613 "clatbs.f"
		} else if (tjj > 0.) {

/*                    0 < abs(A(j,j)) <= SMLNUM: */

#line 617 "clatbs.f"
		    if (xj > tjj * bignum) {

/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
/*                       to avoid overflow when dividing by A(j,j). */

#line 622 "clatbs.f"
			rec = tjj * bignum / xj;
#line 623 "clatbs.f"
			if (cnorm[j] > 1.) {

/*                          Scale by 1/CNORM(j) to avoid overflow when */
/*                          multiplying x(j) times column j. */

#line 628 "clatbs.f"
			    rec /= cnorm[j];
#line 629 "clatbs.f"
			}
#line 630 "clatbs.f"
			csscal_(n, &rec, &x[1], &c__1);
#line 631 "clatbs.f"
			*scale *= rec;
#line 632 "clatbs.f"
			xmax *= rec;
#line 633 "clatbs.f"
		    }
#line 634 "clatbs.f"
		    i__3 = j;
#line 634 "clatbs.f"
		    cladiv_(&z__1, &x[j], &tjjs);
#line 634 "clatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 635 "clatbs.f"
		    i__3 = j;
#line 635 "clatbs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 636 "clatbs.f"
		} else {

/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                    scale = 0, and compute a solution to A*x = 0. */

#line 641 "clatbs.f"
		    i__3 = *n;
#line 641 "clatbs.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 642 "clatbs.f"
			i__4 = i__;
#line 642 "clatbs.f"
			x[i__4].r = 0., x[i__4].i = 0.;
#line 643 "clatbs.f"
/* L100: */
#line 643 "clatbs.f"
		    }
#line 644 "clatbs.f"
		    i__3 = j;
#line 644 "clatbs.f"
		    x[i__3].r = 1., x[i__3].i = 0.;
#line 645 "clatbs.f"
		    xj = 1.;
#line 646 "clatbs.f"
		    *scale = 0.;
#line 647 "clatbs.f"
		    xmax = 0.;
#line 648 "clatbs.f"
		}
#line 649 "clatbs.f"
L105:

/*              Scale x if necessary to avoid overflow when adding a */
/*              multiple of column j of A. */

#line 654 "clatbs.f"
		if (xj > 1.) {
#line 655 "clatbs.f"
		    rec = 1. / xj;
#line 656 "clatbs.f"
		    if (cnorm[j] > (bignum - xmax) * rec) {

/*                    Scale x by 1/(2*abs(x(j))). */

#line 660 "clatbs.f"
			rec *= .5;
#line 661 "clatbs.f"
			csscal_(n, &rec, &x[1], &c__1);
#line 662 "clatbs.f"
			*scale *= rec;
#line 663 "clatbs.f"
		    }
#line 664 "clatbs.f"
		} else if (xj * cnorm[j] > bignum - xmax) {

/*                 Scale x by 1/2. */

#line 668 "clatbs.f"
		    csscal_(n, &c_b36, &x[1], &c__1);
#line 669 "clatbs.f"
		    *scale *= .5;
#line 670 "clatbs.f"
		}

#line 672 "clatbs.f"
		if (upper) {
#line 673 "clatbs.f"
		    if (j > 1) {

/*                    Compute the update */
/*                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) - */
/*                                             x(j)* A(max(1,j-kd):j-1,j) */

/* Computing MIN */
#line 679 "clatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 679 "clatbs.f"
			jlen = min(i__3,i__4);
#line 680 "clatbs.f"
			i__3 = j;
#line 680 "clatbs.f"
			z__2.r = -x[i__3].r, z__2.i = -x[i__3].i;
#line 680 "clatbs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 680 "clatbs.f"
			caxpy_(&jlen, &z__1, &ab[*kd + 1 - jlen + j * ab_dim1]
				, &c__1, &x[j - jlen], &c__1);
#line 682 "clatbs.f"
			i__3 = j - 1;
#line 682 "clatbs.f"
			i__ = icamax_(&i__3, &x[1], &c__1);
#line 683 "clatbs.f"
			i__3 = i__;
#line 683 "clatbs.f"
			xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&x[i__]), abs(d__2));
#line 684 "clatbs.f"
		    }
#line 685 "clatbs.f"
		} else if (j < *n) {

/*                 Compute the update */
/*                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) - */
/*                                          x(j) * A(j+1:min(j+kd,n),j) */

/* Computing MIN */
#line 691 "clatbs.f"
		    i__3 = *kd, i__4 = *n - j;
#line 691 "clatbs.f"
		    jlen = min(i__3,i__4);
#line 692 "clatbs.f"
		    if (jlen > 0) {
#line 692 "clatbs.f"
			i__3 = j;
#line 692 "clatbs.f"
			z__2.r = -x[i__3].r, z__2.i = -x[i__3].i;
#line 692 "clatbs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 692 "clatbs.f"
			caxpy_(&jlen, &z__1, &ab[j * ab_dim1 + 2], &c__1, &x[
				j + 1], &c__1);
#line 692 "clatbs.f"
		    }
#line 695 "clatbs.f"
		    i__3 = *n - j;
#line 695 "clatbs.f"
		    i__ = j + icamax_(&i__3, &x[j + 1], &c__1);
#line 696 "clatbs.f"
		    i__3 = i__;
#line 696 "clatbs.f"
		    xmax = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[
			    i__]), abs(d__2));
#line 697 "clatbs.f"
		}
#line 698 "clatbs.f"
/* L110: */
#line 698 "clatbs.f"
	    }

#line 700 "clatbs.f"
	} else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {

/*           Solve A**T * x = b */

#line 704 "clatbs.f"
	    i__2 = jlast;
#line 704 "clatbs.f"
	    i__1 = jinc;
#line 704 "clatbs.f"
	    for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 709 "clatbs.f"
		i__3 = j;
#line 709 "clatbs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 710 "clatbs.f"
		uscal.r = tscal, uscal.i = 0.;
#line 711 "clatbs.f"
		rec = 1. / max(xmax,1.);
#line 712 "clatbs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 716 "clatbs.f"
		    rec *= .5;
#line 717 "clatbs.f"
		    if (nounit) {
#line 718 "clatbs.f"
			i__3 = maind + j * ab_dim1;
#line 718 "clatbs.f"
			z__1.r = tscal * ab[i__3].r, z__1.i = tscal * ab[i__3]
				.i;
#line 718 "clatbs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 719 "clatbs.f"
		    } else {
#line 720 "clatbs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 721 "clatbs.f"
		    }
#line 722 "clatbs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 723 "clatbs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 727 "clatbs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 727 "clatbs.f"
			rec = min(d__1,d__2);
#line 728 "clatbs.f"
			cladiv_(&z__1, &uscal, &tjjs);
#line 728 "clatbs.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 729 "clatbs.f"
		    }
#line 730 "clatbs.f"
		    if (rec < 1.) {
#line 731 "clatbs.f"
			csscal_(n, &rec, &x[1], &c__1);
#line 732 "clatbs.f"
			*scale *= rec;
#line 733 "clatbs.f"
			xmax *= rec;
#line 734 "clatbs.f"
		    }
#line 735 "clatbs.f"
		}

#line 737 "clatbs.f"
		csumj.r = 0., csumj.i = 0.;
#line 738 "clatbs.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call CDOTU to perform the dot product. */

#line 743 "clatbs.f"
		    if (upper) {
/* Computing MIN */
#line 744 "clatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 744 "clatbs.f"
			jlen = min(i__3,i__4);
#line 745 "clatbs.f"
			cdotu_(&z__1, &jlen, &ab[*kd + 1 - jlen + j * ab_dim1]
				, &c__1, &x[j - jlen], &c__1);
#line 745 "clatbs.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 747 "clatbs.f"
		    } else {
/* Computing MIN */
#line 748 "clatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 748 "clatbs.f"
			jlen = min(i__3,i__4);
#line 749 "clatbs.f"
			if (jlen > 1) {
#line 749 "clatbs.f"
			    cdotu_(&z__1, &jlen, &ab[j * ab_dim1 + 2], &c__1, 
				    &x[j + 1], &c__1);
#line 749 "clatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 749 "clatbs.f"
			}
#line 752 "clatbs.f"
		    }
#line 753 "clatbs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 757 "clatbs.f"
		    if (upper) {
/* Computing MIN */
#line 758 "clatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 758 "clatbs.f"
			jlen = min(i__3,i__4);
#line 759 "clatbs.f"
			i__3 = jlen;
#line 759 "clatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 760 "clatbs.f"
			    i__4 = *kd + i__ - jlen + j * ab_dim1;
#line 760 "clatbs.f"
			    z__3.r = ab[i__4].r * uscal.r - ab[i__4].i * 
				    uscal.i, z__3.i = ab[i__4].r * uscal.i + 
				    ab[i__4].i * uscal.r;
#line 760 "clatbs.f"
			    i__5 = j - jlen - 1 + i__;
#line 760 "clatbs.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 760 "clatbs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 760 "clatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 762 "clatbs.f"
/* L120: */
#line 762 "clatbs.f"
			}
#line 763 "clatbs.f"
		    } else {
/* Computing MIN */
#line 764 "clatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 764 "clatbs.f"
			jlen = min(i__3,i__4);
#line 765 "clatbs.f"
			i__3 = jlen;
#line 765 "clatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 766 "clatbs.f"
			    i__4 = i__ + 1 + j * ab_dim1;
#line 766 "clatbs.f"
			    z__3.r = ab[i__4].r * uscal.r - ab[i__4].i * 
				    uscal.i, z__3.i = ab[i__4].r * uscal.i + 
				    ab[i__4].i * uscal.r;
#line 766 "clatbs.f"
			    i__5 = j + i__;
#line 766 "clatbs.f"
			    z__2.r = z__3.r * x[i__5].r - z__3.i * x[i__5].i, 
				    z__2.i = z__3.r * x[i__5].i + z__3.i * x[
				    i__5].r;
#line 766 "clatbs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 766 "clatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 767 "clatbs.f"
/* L130: */
#line 767 "clatbs.f"
			}
#line 768 "clatbs.f"
		    }
#line 769 "clatbs.f"
		}

#line 771 "clatbs.f"
		z__1.r = tscal, z__1.i = 0.;
#line 771 "clatbs.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 776 "clatbs.f"
		    i__3 = j;
#line 776 "clatbs.f"
		    i__4 = j;
#line 776 "clatbs.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 776 "clatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 777 "clatbs.f"
		    i__3 = j;
#line 777 "clatbs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 778 "clatbs.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 782 "clatbs.f"
			i__3 = maind + j * ab_dim1;
#line 782 "clatbs.f"
			z__1.r = tscal * ab[i__3].r, z__1.i = tscal * ab[i__3]
				.i;
#line 782 "clatbs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 783 "clatbs.f"
		    } else {
#line 784 "clatbs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 785 "clatbs.f"
			if (tscal == 1.) {
#line 785 "clatbs.f"
			    goto L145;
#line 785 "clatbs.f"
			}
#line 787 "clatbs.f"
		    }
#line 788 "clatbs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 789 "clatbs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 793 "clatbs.f"
			if (tjj < 1.) {
#line 794 "clatbs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 798 "clatbs.f"
				rec = 1. / xj;
#line 799 "clatbs.f"
				csscal_(n, &rec, &x[1], &c__1);
#line 800 "clatbs.f"
				*scale *= rec;
#line 801 "clatbs.f"
				xmax *= rec;
#line 802 "clatbs.f"
			    }
#line 803 "clatbs.f"
			}
#line 804 "clatbs.f"
			i__3 = j;
#line 804 "clatbs.f"
			cladiv_(&z__1, &x[j], &tjjs);
#line 804 "clatbs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 805 "clatbs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 809 "clatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 813 "clatbs.f"
			    rec = tjj * bignum / xj;
#line 814 "clatbs.f"
			    csscal_(n, &rec, &x[1], &c__1);
#line 815 "clatbs.f"
			    *scale *= rec;
#line 816 "clatbs.f"
			    xmax *= rec;
#line 817 "clatbs.f"
			}
#line 818 "clatbs.f"
			i__3 = j;
#line 818 "clatbs.f"
			cladiv_(&z__1, &x[j], &tjjs);
#line 818 "clatbs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 819 "clatbs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**T *x = 0. */

#line 824 "clatbs.f"
			i__3 = *n;
#line 824 "clatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 825 "clatbs.f"
			    i__4 = i__;
#line 825 "clatbs.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 826 "clatbs.f"
/* L140: */
#line 826 "clatbs.f"
			}
#line 827 "clatbs.f"
			i__3 = j;
#line 827 "clatbs.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 828 "clatbs.f"
			*scale = 0.;
#line 829 "clatbs.f"
			xmax = 0.;
#line 830 "clatbs.f"
		    }
#line 831 "clatbs.f"
L145:
#line 832 "clatbs.f"
		    ;
#line 832 "clatbs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 837 "clatbs.f"
		    i__3 = j;
#line 837 "clatbs.f"
		    cladiv_(&z__2, &x[j], &tjjs);
#line 837 "clatbs.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 837 "clatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 838 "clatbs.f"
		}
/* Computing MAX */
#line 839 "clatbs.f"
		i__3 = j;
#line 839 "clatbs.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 839 "clatbs.f"
		xmax = max(d__3,d__4);
#line 840 "clatbs.f"
/* L150: */
#line 840 "clatbs.f"
	    }

#line 842 "clatbs.f"
	} else {

/*           Solve A**H * x = b */

#line 846 "clatbs.f"
	    i__1 = jlast;
#line 846 "clatbs.f"
	    i__2 = jinc;
#line 846 "clatbs.f"
	    for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
/*                                    k<>j */

#line 851 "clatbs.f"
		i__3 = j;
#line 851 "clatbs.f"
		xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j]), 
			abs(d__2));
#line 852 "clatbs.f"
		uscal.r = tscal, uscal.i = 0.;
#line 853 "clatbs.f"
		rec = 1. / max(xmax,1.);
#line 854 "clatbs.f"
		if (cnorm[j] > (bignum - xj) * rec) {

/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

#line 858 "clatbs.f"
		    rec *= .5;
#line 859 "clatbs.f"
		    if (nounit) {
#line 860 "clatbs.f"
			d_cnjg(&z__2, &ab[maind + j * ab_dim1]);
#line 860 "clatbs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 860 "clatbs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 861 "clatbs.f"
		    } else {
#line 862 "clatbs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 863 "clatbs.f"
		    }
#line 864 "clatbs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 865 "clatbs.f"
		    if (tjj > 1.) {

/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

/* Computing MIN */
#line 869 "clatbs.f"
			d__1 = 1., d__2 = rec * tjj;
#line 869 "clatbs.f"
			rec = min(d__1,d__2);
#line 870 "clatbs.f"
			cladiv_(&z__1, &uscal, &tjjs);
#line 870 "clatbs.f"
			uscal.r = z__1.r, uscal.i = z__1.i;
#line 871 "clatbs.f"
		    }
#line 872 "clatbs.f"
		    if (rec < 1.) {
#line 873 "clatbs.f"
			csscal_(n, &rec, &x[1], &c__1);
#line 874 "clatbs.f"
			*scale *= rec;
#line 875 "clatbs.f"
			xmax *= rec;
#line 876 "clatbs.f"
		    }
#line 877 "clatbs.f"
		}

#line 879 "clatbs.f"
		csumj.r = 0., csumj.i = 0.;
#line 880 "clatbs.f"
		if (uscal.r == 1. && uscal.i == 0.) {

/*                 If the scaling needed for A in the dot product is 1, */
/*                 call CDOTC to perform the dot product. */

#line 885 "clatbs.f"
		    if (upper) {
/* Computing MIN */
#line 886 "clatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 886 "clatbs.f"
			jlen = min(i__3,i__4);
#line 887 "clatbs.f"
			cdotc_(&z__1, &jlen, &ab[*kd + 1 - jlen + j * ab_dim1]
				, &c__1, &x[j - jlen], &c__1);
#line 887 "clatbs.f"
			csumj.r = z__1.r, csumj.i = z__1.i;
#line 889 "clatbs.f"
		    } else {
/* Computing MIN */
#line 890 "clatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 890 "clatbs.f"
			jlen = min(i__3,i__4);
#line 891 "clatbs.f"
			if (jlen > 1) {
#line 891 "clatbs.f"
			    cdotc_(&z__1, &jlen, &ab[j * ab_dim1 + 2], &c__1, 
				    &x[j + 1], &c__1);
#line 891 "clatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 891 "clatbs.f"
			}
#line 894 "clatbs.f"
		    }
#line 895 "clatbs.f"
		} else {

/*                 Otherwise, use in-line code for the dot product. */

#line 899 "clatbs.f"
		    if (upper) {
/* Computing MIN */
#line 900 "clatbs.f"
			i__3 = *kd, i__4 = j - 1;
#line 900 "clatbs.f"
			jlen = min(i__3,i__4);
#line 901 "clatbs.f"
			i__3 = jlen;
#line 901 "clatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 902 "clatbs.f"
			    d_cnjg(&z__4, &ab[*kd + i__ - jlen + j * ab_dim1])
				    ;
#line 902 "clatbs.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 902 "clatbs.f"
			    i__4 = j - jlen - 1 + i__;
#line 902 "clatbs.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 902 "clatbs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 902 "clatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 904 "clatbs.f"
/* L160: */
#line 904 "clatbs.f"
			}
#line 905 "clatbs.f"
		    } else {
/* Computing MIN */
#line 906 "clatbs.f"
			i__3 = *kd, i__4 = *n - j;
#line 906 "clatbs.f"
			jlen = min(i__3,i__4);
#line 907 "clatbs.f"
			i__3 = jlen;
#line 907 "clatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 908 "clatbs.f"
			    d_cnjg(&z__4, &ab[i__ + 1 + j * ab_dim1]);
#line 908 "clatbs.f"
			    z__3.r = z__4.r * uscal.r - z__4.i * uscal.i, 
				    z__3.i = z__4.r * uscal.i + z__4.i * 
				    uscal.r;
#line 908 "clatbs.f"
			    i__4 = j + i__;
#line 908 "clatbs.f"
			    z__2.r = z__3.r * x[i__4].r - z__3.i * x[i__4].i, 
				    z__2.i = z__3.r * x[i__4].i + z__3.i * x[
				    i__4].r;
#line 908 "clatbs.f"
			    z__1.r = csumj.r + z__2.r, z__1.i = csumj.i + 
				    z__2.i;
#line 908 "clatbs.f"
			    csumj.r = z__1.r, csumj.i = z__1.i;
#line 910 "clatbs.f"
/* L170: */
#line 910 "clatbs.f"
			}
#line 911 "clatbs.f"
		    }
#line 912 "clatbs.f"
		}

#line 914 "clatbs.f"
		z__1.r = tscal, z__1.i = 0.;
#line 914 "clatbs.f"
		if (uscal.r == z__1.r && uscal.i == z__1.i) {

/*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j) */
/*                 was not used to scale the dotproduct. */

#line 919 "clatbs.f"
		    i__3 = j;
#line 919 "clatbs.f"
		    i__4 = j;
#line 919 "clatbs.f"
		    z__1.r = x[i__4].r - csumj.r, z__1.i = x[i__4].i - 
			    csumj.i;
#line 919 "clatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 920 "clatbs.f"
		    i__3 = j;
#line 920 "clatbs.f"
		    xj = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[j])
			    , abs(d__2));
#line 921 "clatbs.f"
		    if (nounit) {

/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

#line 925 "clatbs.f"
			d_cnjg(&z__2, &ab[maind + j * ab_dim1]);
#line 925 "clatbs.f"
			z__1.r = tscal * z__2.r, z__1.i = tscal * z__2.i;
#line 925 "clatbs.f"
			tjjs.r = z__1.r, tjjs.i = z__1.i;
#line 926 "clatbs.f"
		    } else {
#line 927 "clatbs.f"
			tjjs.r = tscal, tjjs.i = 0.;
#line 928 "clatbs.f"
			if (tscal == 1.) {
#line 928 "clatbs.f"
			    goto L185;
#line 928 "clatbs.f"
			}
#line 930 "clatbs.f"
		    }
#line 931 "clatbs.f"
		    tjj = (d__1 = tjjs.r, abs(d__1)) + (d__2 = d_imag(&tjjs), 
			    abs(d__2));
#line 932 "clatbs.f"
		    if (tjj > smlnum) {

/*                       abs(A(j,j)) > SMLNUM: */

#line 936 "clatbs.f"
			if (tjj < 1.) {
#line 937 "clatbs.f"
			    if (xj > tjj * bignum) {

/*                             Scale X by 1/abs(x(j)). */

#line 941 "clatbs.f"
				rec = 1. / xj;
#line 942 "clatbs.f"
				csscal_(n, &rec, &x[1], &c__1);
#line 943 "clatbs.f"
				*scale *= rec;
#line 944 "clatbs.f"
				xmax *= rec;
#line 945 "clatbs.f"
			    }
#line 946 "clatbs.f"
			}
#line 947 "clatbs.f"
			i__3 = j;
#line 947 "clatbs.f"
			cladiv_(&z__1, &x[j], &tjjs);
#line 947 "clatbs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 948 "clatbs.f"
		    } else if (tjj > 0.) {

/*                       0 < abs(A(j,j)) <= SMLNUM: */

#line 952 "clatbs.f"
			if (xj > tjj * bignum) {

/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

#line 956 "clatbs.f"
			    rec = tjj * bignum / xj;
#line 957 "clatbs.f"
			    csscal_(n, &rec, &x[1], &c__1);
#line 958 "clatbs.f"
			    *scale *= rec;
#line 959 "clatbs.f"
			    xmax *= rec;
#line 960 "clatbs.f"
			}
#line 961 "clatbs.f"
			i__3 = j;
#line 961 "clatbs.f"
			cladiv_(&z__1, &x[j], &tjjs);
#line 961 "clatbs.f"
			x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 962 "clatbs.f"
		    } else {

/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
/*                       scale = 0 and compute a solution to A**H *x = 0. */

#line 967 "clatbs.f"
			i__3 = *n;
#line 967 "clatbs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 968 "clatbs.f"
			    i__4 = i__;
#line 968 "clatbs.f"
			    x[i__4].r = 0., x[i__4].i = 0.;
#line 969 "clatbs.f"
/* L180: */
#line 969 "clatbs.f"
			}
#line 970 "clatbs.f"
			i__3 = j;
#line 970 "clatbs.f"
			x[i__3].r = 1., x[i__3].i = 0.;
#line 971 "clatbs.f"
			*scale = 0.;
#line 972 "clatbs.f"
			xmax = 0.;
#line 973 "clatbs.f"
		    }
#line 974 "clatbs.f"
L185:
#line 975 "clatbs.f"
		    ;
#line 975 "clatbs.f"
		} else {

/*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot */
/*                 product has already been divided by 1/A(j,j). */

#line 980 "clatbs.f"
		    i__3 = j;
#line 980 "clatbs.f"
		    cladiv_(&z__2, &x[j], &tjjs);
#line 980 "clatbs.f"
		    z__1.r = z__2.r - csumj.r, z__1.i = z__2.i - csumj.i;
#line 980 "clatbs.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 981 "clatbs.f"
		}
/* Computing MAX */
#line 982 "clatbs.f"
		i__3 = j;
#line 982 "clatbs.f"
		d__3 = xmax, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&x[j]), abs(d__2));
#line 982 "clatbs.f"
		xmax = max(d__3,d__4);
#line 983 "clatbs.f"
/* L190: */
#line 983 "clatbs.f"
	    }
#line 984 "clatbs.f"
	}
#line 985 "clatbs.f"
	*scale /= tscal;
#line 986 "clatbs.f"
    }

/*     Scale the column norms by 1/TSCAL for return. */

#line 990 "clatbs.f"
    if (tscal != 1.) {
#line 991 "clatbs.f"
	d__1 = 1. / tscal;
#line 991 "clatbs.f"
	sscal_(n, &d__1, &cnorm[1], &c__1);
#line 992 "clatbs.f"
    }

#line 994 "clatbs.f"
    return 0;

/*     End of CLATBS */

} /* clatbs_ */

