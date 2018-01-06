#line 1 "dsytf2_rook.f"
/* dsytf2_rook.f -- translated by f2c (version 20100827).
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

#line 1 "dsytf2_rook.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DSYTF2_ROOK computes the factorization of a real symmetric indefinite matrix using the bounded 
Bunch-Kaufman ("rook") diagonal pivoting method (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTF2_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytf2_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytf2_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytf2_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTF2_ROOK( UPLO, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTF2_ROOK computes the factorization of a real symmetric matrix A */
/* > using the bounded Bunch-Kaufman ("rook") diagonal pivoting method: */
/* > */
/* >    A = U*D*U**T  or  A = L*D*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, U**T is the transpose of U, and D is symmetric and */
/* > block diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored: */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          n-by-n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n-by-n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, the block diagonal matrix D and the multipliers used */
/* >          to obtain the factor U or L (see below for further details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D. */
/* > */
/* >          If UPLO = 'U': */
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) */
/* >             were interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and */
/* >             columns k and -IPIV(k) were interchanged and rows and */
/* >             columns k-1 and -IPIV(k-1) were inerchaged, */
/* >             D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
/* > */
/* >          If UPLO = 'L': */
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) */
/* >             were interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and */
/* >             columns k and -IPIV(k) were interchanged and rows and */
/* >             columns k+1 and -IPIV(k+1) were inerchaged, */
/* >             D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* >          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization */
/* >               has been completed, but the block diagonal matrix D is */
/* >               exactly singular, and division by zero will occur if it */
/* >               is used to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2013 */

/* > \ingroup doubleSYcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  If UPLO = 'U', then A = U*D*U**T, where */
/* >     U = P(n)*U(n)* ... *P(k)U(k)* ..., */
/* >  i.e., U is a product of terms P(k)*U(k), where k decreases from n to */
/* >  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1 */
/* >  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as */
/* >  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such */
/* >  that if the diagonal block D(k) is of order s (s = 1 or 2), then */
/* > */
/* >             (   I    v    0   )   k-s */
/* >     U(k) =  (   0    I    0   )   s */
/* >             (   0    0    I   )   n-k */
/* >                k-s   s   n-k */
/* > */
/* >  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k). */
/* >  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k), */
/* >  and A(k,k), and v overwrites A(1:k-2,k-1:k). */
/* > */
/* >  If UPLO = 'L', then A = L*D*L**T, where */
/* >     L = P(1)*L(1)* ... *P(k)*L(k)* ..., */
/* >  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to */
/* >  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1 */
/* >  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as */
/* >  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such */
/* >  that if the diagonal block D(k) is of order s (s = 1 or 2), then */
/* > */
/* >             (   I    0     0   )  k-1 */
/* >     L(k) =  (   0    I     0   )  s */
/* >             (   0    v     I   )  n-k-s+1 */
/* >                k-1   s  n-k-s+1 */
/* > */
/* >  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k). */
/* >  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k), */
/* >  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1). */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2013,     Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* >  01-01-96 - Based on modifications by */
/* >    J. Lewis, Boeing Computer Services Company */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville abd , USA */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int dsytf2_rook__(char *uplo, integer *n, doublereal *a, 
	integer *lda, integer *ipiv, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, p;
    static doublereal t, d11, d12, d21, d22;
    static integer ii, kk, kp;
    static doublereal wk, wkm1, wkp1;
    static logical done;
    static integer imax, jmax;
    extern /* Subroutine */ int dsyr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dtemp, sfmin;
    static integer itemp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer kstep;
    static logical upper;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal absakk;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal colmax, rowmax;


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

/*     Test the input parameters. */

#line 242 "dsytf2_rook.f"
    /* Parameter adjustments */
#line 242 "dsytf2_rook.f"
    a_dim1 = *lda;
#line 242 "dsytf2_rook.f"
    a_offset = 1 + a_dim1;
#line 242 "dsytf2_rook.f"
    a -= a_offset;
#line 242 "dsytf2_rook.f"
    --ipiv;
#line 242 "dsytf2_rook.f"

#line 242 "dsytf2_rook.f"
    /* Function Body */
#line 242 "dsytf2_rook.f"
    *info = 0;
#line 243 "dsytf2_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 244 "dsytf2_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 245 "dsytf2_rook.f"
	*info = -1;
#line 246 "dsytf2_rook.f"
    } else if (*n < 0) {
#line 247 "dsytf2_rook.f"
	*info = -2;
#line 248 "dsytf2_rook.f"
    } else if (*lda < max(1,*n)) {
#line 249 "dsytf2_rook.f"
	*info = -4;
#line 250 "dsytf2_rook.f"
    }
#line 251 "dsytf2_rook.f"
    if (*info != 0) {
#line 252 "dsytf2_rook.f"
	i__1 = -(*info);
#line 252 "dsytf2_rook.f"
	xerbla_("DSYTF2_ROOK", &i__1, (ftnlen)11);
#line 253 "dsytf2_rook.f"
	return 0;
#line 254 "dsytf2_rook.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 258 "dsytf2_rook.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 262 "dsytf2_rook.f"
    sfmin = dlamch_("S", (ftnlen)1);

#line 264 "dsytf2_rook.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 271 "dsytf2_rook.f"
	k = *n;
#line 272 "dsytf2_rook.f"
L10:

/*        If K < 1, exit from loop */

#line 276 "dsytf2_rook.f"
	if (k < 1) {
#line 276 "dsytf2_rook.f"
	    goto L70;
#line 276 "dsytf2_rook.f"
	}
#line 278 "dsytf2_rook.f"
	kstep = 1;
#line 279 "dsytf2_rook.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 284 "dsytf2_rook.f"
	absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 290 "dsytf2_rook.f"
	if (k > 1) {
#line 291 "dsytf2_rook.f"
	    i__1 = k - 1;
#line 291 "dsytf2_rook.f"
	    imax = idamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 292 "dsytf2_rook.f"
	    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
#line 293 "dsytf2_rook.f"
	} else {
#line 294 "dsytf2_rook.f"
	    colmax = 0.;
#line 295 "dsytf2_rook.f"
	}

#line 297 "dsytf2_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 301 "dsytf2_rook.f"
	    if (*info == 0) {
#line 301 "dsytf2_rook.f"
		*info = k;
#line 301 "dsytf2_rook.f"
	    }
#line 303 "dsytf2_rook.f"
	    kp = k;
#line 304 "dsytf2_rook.f"
	} else {

/*           Test for interchange */

/*           Equivalent to testing for (used to handle NaN and Inf) */
/*           ABSAKK.GE.ALPHA*COLMAX */

#line 311 "dsytf2_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, */
/*              use 1-by-1 pivot block */

#line 316 "dsytf2_rook.f"
		kp = k;
#line 317 "dsytf2_rook.f"
	    } else {

#line 319 "dsytf2_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 323 "dsytf2_rook.f"
L12:

/*                 Begin pivot search loop body */

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 331 "dsytf2_rook.f"
		if (imax != k) {
#line 332 "dsytf2_rook.f"
		    i__1 = k - imax;
#line 332 "dsytf2_rook.f"
		    jmax = imax + idamax_(&i__1, &a[imax + (imax + 1) * 
			    a_dim1], lda);
#line 334 "dsytf2_rook.f"
		    rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
#line 335 "dsytf2_rook.f"
		} else {
#line 336 "dsytf2_rook.f"
		    rowmax = 0.;
#line 337 "dsytf2_rook.f"
		}

#line 339 "dsytf2_rook.f"
		if (imax > 1) {
#line 340 "dsytf2_rook.f"
		    i__1 = imax - 1;
#line 340 "dsytf2_rook.f"
		    itemp = idamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
#line 341 "dsytf2_rook.f"
		    dtemp = (d__1 = a[itemp + imax * a_dim1], abs(d__1));
#line 342 "dsytf2_rook.f"
		    if (dtemp > rowmax) {
#line 343 "dsytf2_rook.f"
			rowmax = dtemp;
#line 344 "dsytf2_rook.f"
			jmax = itemp;
#line 345 "dsytf2_rook.f"
		    }
#line 346 "dsytf2_rook.f"
		}

/*                 Equivalent to testing for (used to handle NaN and Inf) */
/*                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */

#line 351 "dsytf2_rook.f"
		if (! ((d__1 = a[imax + imax * a_dim1], abs(d__1)) < alpha * 
			rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 357 "dsytf2_rook.f"
		    kp = imax;
#line 358 "dsytf2_rook.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX .EQ. COLMAX, */
/*                 used to handle NaN and Inf */

#line 363 "dsytf2_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 368 "dsytf2_rook.f"
		    kp = imax;
#line 369 "dsytf2_rook.f"
		    kstep = 2;
#line 370 "dsytf2_rook.f"
		    done = TRUE_;
#line 371 "dsytf2_rook.f"
		} else {

/*                    Pivot NOT found, set variables and repeat */

#line 375 "dsytf2_rook.f"
		    p = imax;
#line 376 "dsytf2_rook.f"
		    colmax = rowmax;
#line 377 "dsytf2_rook.f"
		    imax = jmax;
#line 378 "dsytf2_rook.f"
		}

/*                 End pivot search loop body */

#line 382 "dsytf2_rook.f"
		if (! done) {
#line 382 "dsytf2_rook.f"
		    goto L12;
#line 382 "dsytf2_rook.f"
		}

#line 384 "dsytf2_rook.f"
	    }

/*           Swap TWO rows and TWO columns */

/*           First swap */

#line 390 "dsytf2_rook.f"
	    if (kstep == 2 && p != k) {

/*              Interchange rows and column K and P in the leading */
/*              submatrix A(1:k,1:k) if we have a 2-by-2 pivot */

#line 395 "dsytf2_rook.f"
		if (p > 1) {
#line 395 "dsytf2_rook.f"
		    i__1 = p - 1;
#line 395 "dsytf2_rook.f"
		    dswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 
			    1], &c__1);
#line 395 "dsytf2_rook.f"
		}
#line 397 "dsytf2_rook.f"
		if (p < k - 1) {
#line 397 "dsytf2_rook.f"
		    i__1 = k - p - 1;
#line 397 "dsytf2_rook.f"
		    dswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 
			    1) * a_dim1], lda);
#line 397 "dsytf2_rook.f"
		}
#line 400 "dsytf2_rook.f"
		t = a[k + k * a_dim1];
#line 401 "dsytf2_rook.f"
		a[k + k * a_dim1] = a[p + p * a_dim1];
#line 402 "dsytf2_rook.f"
		a[p + p * a_dim1] = t;
#line 403 "dsytf2_rook.f"
	    }

/*           Second swap */

#line 407 "dsytf2_rook.f"
	    kk = k - kstep + 1;
#line 408 "dsytf2_rook.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 413 "dsytf2_rook.f"
		if (kp > 1) {
#line 413 "dsytf2_rook.f"
		    i__1 = kp - 1;
#line 413 "dsytf2_rook.f"
		    dswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 413 "dsytf2_rook.f"
		}
#line 415 "dsytf2_rook.f"
		if (kk > 1 && kp < kk - 1) {
#line 415 "dsytf2_rook.f"
		    i__1 = kk - kp - 1;
#line 415 "dsytf2_rook.f"
		    dswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (
			    kp + 1) * a_dim1], lda);
#line 415 "dsytf2_rook.f"
		}
#line 418 "dsytf2_rook.f"
		t = a[kk + kk * a_dim1];
#line 419 "dsytf2_rook.f"
		a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
#line 420 "dsytf2_rook.f"
		a[kp + kp * a_dim1] = t;
#line 421 "dsytf2_rook.f"
		if (kstep == 2) {
#line 422 "dsytf2_rook.f"
		    t = a[k - 1 + k * a_dim1];
#line 423 "dsytf2_rook.f"
		    a[k - 1 + k * a_dim1] = a[kp + k * a_dim1];
#line 424 "dsytf2_rook.f"
		    a[kp + k * a_dim1] = t;
#line 425 "dsytf2_rook.f"
		}
#line 426 "dsytf2_rook.f"
	    }

/*           Update the leading submatrix */

#line 430 "dsytf2_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

#line 438 "dsytf2_rook.f"
		if (k > 1) {

/*                 Perform a rank-1 update of A(1:k-1,1:k-1) and */
/*                 store U(k) in column k */

#line 443 "dsytf2_rook.f"
		    if ((d__1 = a[k + k * a_dim1], abs(d__1)) >= sfmin) {

/*                    Perform a rank-1 update of A(1:k-1,1:k-1) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*1/D(k)*W(k)**T */

#line 449 "dsytf2_rook.f"
			d11 = 1. / a[k + k * a_dim1];
#line 450 "dsytf2_rook.f"
			i__1 = k - 1;
#line 450 "dsytf2_rook.f"
			d__1 = -d11;
#line 450 "dsytf2_rook.f"
			dsyr_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);

/*                    Store U(k) in column k */

#line 454 "dsytf2_rook.f"
			i__1 = k - 1;
#line 454 "dsytf2_rook.f"
			dscal_(&i__1, &d11, &a[k * a_dim1 + 1], &c__1);
#line 455 "dsytf2_rook.f"
		    } else {

/*                    Store L(k) in column K */

#line 459 "dsytf2_rook.f"
			d11 = a[k + k * a_dim1];
#line 460 "dsytf2_rook.f"
			i__1 = k - 1;
#line 460 "dsytf2_rook.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 461 "dsytf2_rook.f"
			    a[ii + k * a_dim1] /= d11;
#line 462 "dsytf2_rook.f"
/* L16: */
#line 462 "dsytf2_rook.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 469 "dsytf2_rook.f"
			i__1 = k - 1;
#line 469 "dsytf2_rook.f"
			d__1 = -d11;
#line 469 "dsytf2_rook.f"
			dsyr_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);
#line 470 "dsytf2_rook.f"
		    }
#line 471 "dsytf2_rook.f"
		}

#line 473 "dsytf2_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 489 "dsytf2_rook.f"
		if (k > 2) {

#line 491 "dsytf2_rook.f"
		    d12 = a[k - 1 + k * a_dim1];
#line 492 "dsytf2_rook.f"
		    d22 = a[k - 1 + (k - 1) * a_dim1] / d12;
#line 493 "dsytf2_rook.f"
		    d11 = a[k + k * a_dim1] / d12;
#line 494 "dsytf2_rook.f"
		    t = 1. / (d11 * d22 - 1.);

#line 496 "dsytf2_rook.f"
		    for (j = k - 2; j >= 1; --j) {

#line 498 "dsytf2_rook.f"
			wkm1 = t * (d11 * a[j + (k - 1) * a_dim1] - a[j + k * 
				a_dim1]);
#line 499 "dsytf2_rook.f"
			wk = t * (d22 * a[j + k * a_dim1] - a[j + (k - 1) * 
				a_dim1]);

#line 501 "dsytf2_rook.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 502 "dsytf2_rook.f"
			    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ 
				    + k * a_dim1] / d12 * wk - a[i__ + (k - 1)
				     * a_dim1] / d12 * wkm1;
#line 504 "dsytf2_rook.f"
/* L20: */
#line 504 "dsytf2_rook.f"
			}

/*                    Store U(k) and U(k-1) in cols k and k-1 for row J */

#line 508 "dsytf2_rook.f"
			a[j + k * a_dim1] = wk / d12;
#line 509 "dsytf2_rook.f"
			a[j + (k - 1) * a_dim1] = wkm1 / d12;

#line 511 "dsytf2_rook.f"
/* L30: */
#line 511 "dsytf2_rook.f"
		    }

#line 513 "dsytf2_rook.f"
		}

#line 515 "dsytf2_rook.f"
	    }
#line 516 "dsytf2_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 520 "dsytf2_rook.f"
	if (kstep == 1) {
#line 521 "dsytf2_rook.f"
	    ipiv[k] = kp;
#line 522 "dsytf2_rook.f"
	} else {
#line 523 "dsytf2_rook.f"
	    ipiv[k] = -p;
#line 524 "dsytf2_rook.f"
	    ipiv[k - 1] = -kp;
#line 525 "dsytf2_rook.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 529 "dsytf2_rook.f"
	k -= kstep;
#line 530 "dsytf2_rook.f"
	goto L10;

#line 532 "dsytf2_rook.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 539 "dsytf2_rook.f"
	k = 1;
#line 540 "dsytf2_rook.f"
L40:

/*        If K > N, exit from loop */

#line 544 "dsytf2_rook.f"
	if (k > *n) {
#line 544 "dsytf2_rook.f"
	    goto L70;
#line 544 "dsytf2_rook.f"
	}
#line 546 "dsytf2_rook.f"
	kstep = 1;
#line 547 "dsytf2_rook.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 552 "dsytf2_rook.f"
	absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 558 "dsytf2_rook.f"
	if (k < *n) {
#line 559 "dsytf2_rook.f"
	    i__1 = *n - k;
#line 559 "dsytf2_rook.f"
	    imax = k + idamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 560 "dsytf2_rook.f"
	    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
#line 561 "dsytf2_rook.f"
	} else {
#line 562 "dsytf2_rook.f"
	    colmax = 0.;
#line 563 "dsytf2_rook.f"
	}

#line 565 "dsytf2_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 569 "dsytf2_rook.f"
	    if (*info == 0) {
#line 569 "dsytf2_rook.f"
		*info = k;
#line 569 "dsytf2_rook.f"
	    }
#line 571 "dsytf2_rook.f"
	    kp = k;
#line 572 "dsytf2_rook.f"
	} else {

/*           Test for interchange */

/*           Equivalent to testing for (used to handle NaN and Inf) */
/*           ABSAKK.GE.ALPHA*COLMAX */

#line 579 "dsytf2_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 583 "dsytf2_rook.f"
		kp = k;
#line 584 "dsytf2_rook.f"
	    } else {

#line 586 "dsytf2_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 590 "dsytf2_rook.f"
L42:

/*                 Begin pivot search loop body */

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 598 "dsytf2_rook.f"
		if (imax != k) {
#line 599 "dsytf2_rook.f"
		    i__1 = imax - k;
#line 599 "dsytf2_rook.f"
		    jmax = k - 1 + idamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 600 "dsytf2_rook.f"
		    rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
#line 601 "dsytf2_rook.f"
		} else {
#line 602 "dsytf2_rook.f"
		    rowmax = 0.;
#line 603 "dsytf2_rook.f"
		}

#line 605 "dsytf2_rook.f"
		if (imax < *n) {
#line 606 "dsytf2_rook.f"
		    i__1 = *n - imax;
#line 606 "dsytf2_rook.f"
		    itemp = imax + idamax_(&i__1, &a[imax + 1 + imax * a_dim1]
			    , &c__1);
#line 608 "dsytf2_rook.f"
		    dtemp = (d__1 = a[itemp + imax * a_dim1], abs(d__1));
#line 609 "dsytf2_rook.f"
		    if (dtemp > rowmax) {
#line 610 "dsytf2_rook.f"
			rowmax = dtemp;
#line 611 "dsytf2_rook.f"
			jmax = itemp;
#line 612 "dsytf2_rook.f"
		    }
#line 613 "dsytf2_rook.f"
		}

/*                 Equivalent to testing for (used to handle NaN and Inf) */
/*                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */

#line 618 "dsytf2_rook.f"
		if (! ((d__1 = a[imax + imax * a_dim1], abs(d__1)) < alpha * 
			rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 624 "dsytf2_rook.f"
		    kp = imax;
#line 625 "dsytf2_rook.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX .EQ. COLMAX, */
/*                 used to handle NaN and Inf */

#line 630 "dsytf2_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 635 "dsytf2_rook.f"
		    kp = imax;
#line 636 "dsytf2_rook.f"
		    kstep = 2;
#line 637 "dsytf2_rook.f"
		    done = TRUE_;
#line 638 "dsytf2_rook.f"
		} else {

/*                    Pivot NOT found, set variables and repeat */

#line 642 "dsytf2_rook.f"
		    p = imax;
#line 643 "dsytf2_rook.f"
		    colmax = rowmax;
#line 644 "dsytf2_rook.f"
		    imax = jmax;
#line 645 "dsytf2_rook.f"
		}

/*                 End pivot search loop body */

#line 649 "dsytf2_rook.f"
		if (! done) {
#line 649 "dsytf2_rook.f"
		    goto L42;
#line 649 "dsytf2_rook.f"
		}

#line 651 "dsytf2_rook.f"
	    }

/*           Swap TWO rows and TWO columns */

/*           First swap */

#line 657 "dsytf2_rook.f"
	    if (kstep == 2 && p != k) {

/*              Interchange rows and column K and P in the trailing */
/*              submatrix A(k:n,k:n) if we have a 2-by-2 pivot */

#line 662 "dsytf2_rook.f"
		if (p < *n) {
#line 662 "dsytf2_rook.f"
		    i__1 = *n - p;
#line 662 "dsytf2_rook.f"
		    dswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p 
			    * a_dim1], &c__1);
#line 662 "dsytf2_rook.f"
		}
#line 664 "dsytf2_rook.f"
		if (p > k + 1) {
#line 664 "dsytf2_rook.f"
		    i__1 = p - k - 1;
#line 664 "dsytf2_rook.f"
		    dswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[p + (k + 
			    1) * a_dim1], lda);
#line 664 "dsytf2_rook.f"
		}
#line 666 "dsytf2_rook.f"
		t = a[k + k * a_dim1];
#line 667 "dsytf2_rook.f"
		a[k + k * a_dim1] = a[p + p * a_dim1];
#line 668 "dsytf2_rook.f"
		a[p + p * a_dim1] = t;
#line 669 "dsytf2_rook.f"
	    }

/*           Second swap */

#line 673 "dsytf2_rook.f"
	    kk = k + kstep - 1;
#line 674 "dsytf2_rook.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 679 "dsytf2_rook.f"
		if (kp < *n) {
#line 679 "dsytf2_rook.f"
		    i__1 = *n - kp;
#line 679 "dsytf2_rook.f"
		    dswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 679 "dsytf2_rook.f"
		}
#line 681 "dsytf2_rook.f"
		if (kk < *n && kp > kk + 1) {
#line 681 "dsytf2_rook.f"
		    i__1 = kp - kk - 1;
#line 681 "dsytf2_rook.f"
		    dswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (
			    kk + 1) * a_dim1], lda);
#line 681 "dsytf2_rook.f"
		}
#line 684 "dsytf2_rook.f"
		t = a[kk + kk * a_dim1];
#line 685 "dsytf2_rook.f"
		a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
#line 686 "dsytf2_rook.f"
		a[kp + kp * a_dim1] = t;
#line 687 "dsytf2_rook.f"
		if (kstep == 2) {
#line 688 "dsytf2_rook.f"
		    t = a[k + 1 + k * a_dim1];
#line 689 "dsytf2_rook.f"
		    a[k + 1 + k * a_dim1] = a[kp + k * a_dim1];
#line 690 "dsytf2_rook.f"
		    a[kp + k * a_dim1] = t;
#line 691 "dsytf2_rook.f"
		}
#line 692 "dsytf2_rook.f"
	    }

/*           Update the trailing submatrix */

#line 696 "dsytf2_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 704 "dsytf2_rook.f"
		if (k < *n) {

/*              Perform a rank-1 update of A(k+1:n,k+1:n) and */
/*              store L(k) in column k */

#line 709 "dsytf2_rook.f"
		    if ((d__1 = a[k + k * a_dim1], abs(d__1)) >= sfmin) {

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */

#line 715 "dsytf2_rook.f"
			d11 = 1. / a[k + k * a_dim1];
#line 716 "dsytf2_rook.f"
			i__1 = *n - k;
#line 716 "dsytf2_rook.f"
			d__1 = -d11;
#line 716 "dsytf2_rook.f"
			dsyr_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);

/*                    Store L(k) in column k */

#line 721 "dsytf2_rook.f"
			i__1 = *n - k;
#line 721 "dsytf2_rook.f"
			dscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
#line 722 "dsytf2_rook.f"
		    } else {

/*                    Store L(k) in column k */

#line 726 "dsytf2_rook.f"
			d11 = a[k + k * a_dim1];
#line 727 "dsytf2_rook.f"
			i__1 = *n;
#line 727 "dsytf2_rook.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 728 "dsytf2_rook.f"
			    a[ii + k * a_dim1] /= d11;
#line 729 "dsytf2_rook.f"
/* L46: */
#line 729 "dsytf2_rook.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 736 "dsytf2_rook.f"
			i__1 = *n - k;
#line 736 "dsytf2_rook.f"
			d__1 = -d11;
#line 736 "dsytf2_rook.f"
			dsyr_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);
#line 738 "dsytf2_rook.f"
		    }
#line 739 "dsytf2_rook.f"
		}

#line 741 "dsytf2_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */


/*              Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T */
/*                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 758 "dsytf2_rook.f"
		if (k < *n - 1) {

#line 760 "dsytf2_rook.f"
		    d21 = a[k + 1 + k * a_dim1];
#line 761 "dsytf2_rook.f"
		    d11 = a[k + 1 + (k + 1) * a_dim1] / d21;
#line 762 "dsytf2_rook.f"
		    d22 = a[k + k * a_dim1] / d21;
#line 763 "dsytf2_rook.f"
		    t = 1. / (d11 * d22 - 1.);

#line 765 "dsytf2_rook.f"
		    i__1 = *n;
#line 765 "dsytf2_rook.f"
		    for (j = k + 2; j <= i__1; ++j) {

/*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */

#line 769 "dsytf2_rook.f"
			wk = t * (d11 * a[j + k * a_dim1] - a[j + (k + 1) * 
				a_dim1]);
#line 770 "dsytf2_rook.f"
			wkp1 = t * (d22 * a[j + (k + 1) * a_dim1] - a[j + k * 
				a_dim1]);

/*                    Perform a rank-2 update of A(k+2:n,k+2:n) */

#line 774 "dsytf2_rook.f"
			i__2 = *n;
#line 774 "dsytf2_rook.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 775 "dsytf2_rook.f"
			    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ 
				    + k * a_dim1] / d21 * wk - a[i__ + (k + 1)
				     * a_dim1] / d21 * wkp1;
#line 777 "dsytf2_rook.f"
/* L50: */
#line 777 "dsytf2_rook.f"
			}

/*                    Store L(k) and L(k+1) in cols k and k+1 for row J */

#line 781 "dsytf2_rook.f"
			a[j + k * a_dim1] = wk / d21;
#line 782 "dsytf2_rook.f"
			a[j + (k + 1) * a_dim1] = wkp1 / d21;

#line 784 "dsytf2_rook.f"
/* L60: */
#line 784 "dsytf2_rook.f"
		    }

#line 786 "dsytf2_rook.f"
		}

#line 788 "dsytf2_rook.f"
	    }
#line 789 "dsytf2_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 793 "dsytf2_rook.f"
	if (kstep == 1) {
#line 794 "dsytf2_rook.f"
	    ipiv[k] = kp;
#line 795 "dsytf2_rook.f"
	} else {
#line 796 "dsytf2_rook.f"
	    ipiv[k] = -p;
#line 797 "dsytf2_rook.f"
	    ipiv[k + 1] = -kp;
#line 798 "dsytf2_rook.f"
	}

/*        Increase K and return to the start of the main loop */

#line 802 "dsytf2_rook.f"
	k += kstep;
#line 803 "dsytf2_rook.f"
	goto L40;

#line 805 "dsytf2_rook.f"
    }

#line 807 "dsytf2_rook.f"
L70:

#line 809 "dsytf2_rook.f"
    return 0;

/*     End of DSYTF2_ROOK */

} /* dsytf2_rook__ */

