#line 1 "dsytf2_rk.f"
/* dsytf2_rk.f -- translated by f2c (version 20100827).
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

#line 1 "dsytf2_rk.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DSYTF2_RK computes the factorization of a real symmetric indefinite matrix using the bounded Bu
nch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTF2_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytf2_
rk.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytf2_
rk.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytf2_
rk.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), E ( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > DSYTF2_RK computes the factorization of a real symmetric matrix A */
/* > using the bounded Bunch-Kaufman (rook) diagonal pivoting method: */
/* > */
/* >    A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is symmetric and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
/* > For more information see Further Details section. */
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
/* >          On entry, the symmetric matrix A. */
/* >            If UPLO = 'U': the leading N-by-N upper triangular part */
/* >            of A contains the upper triangular part of the matrix A, */
/* >            and the strictly lower triangular part of A is not */
/* >            referenced. */
/* > */
/* >            If UPLO = 'L': the leading N-by-N lower triangular part */
/* >            of A contains the lower triangular part of the matrix A, */
/* >            and the strictly upper triangular part of A is not */
/* >            referenced. */
/* > */
/* >          On exit, contains: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N) */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          NOTE: For 1-by-1 diagonal block D(k), where */
/* >          1 <= k <= N, the element E(k) is set to 0 in both */
/* >          UPLO = 'U' or UPLO = 'L' cases. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          IPIV describes the permutation matrix P in the factorization */
/* >          of matrix A as follows. The absolute value of IPIV(k) */
/* >          represents the index of row and column that were */
/* >          interchanged with the k-th row and column. The value of UPLO */
/* >          describes the order in which the interchanges were applied. */
/* >          Also, the sign of IPIV represents the block structure of */
/* >          the symmetric block diagonal matrix D with 1-by-1 or 2-by-2 */
/* >          diagonal blocks which correspond to 1 or 2 interchanges */
/* >          at each factorization step. For more info see Further */
/* >          Details section. */
/* > */
/* >          If UPLO = 'U', */
/* >          ( in factorization order, k decreases from N to 1 ): */
/* >            a) A single positive entry IPIV(k) > 0 means: */
/* >               D(k,k) is a 1-by-1 diagonal block. */
/* >               If IPIV(k) != k, rows and columns k and IPIV(k) were */
/* >               interchanged in the matrix A(1:N,1:N); */
/* >               If IPIV(k) = k, no interchange occurred. */
/* > */
/* >            b) A pair of consecutive negative entries */
/* >               IPIV(k) < 0 and IPIV(k-1) < 0 means: */
/* >               D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
/* >               (NOTE: negative entries in IPIV appear ONLY in pairs). */
/* >               1) If -IPIV(k) != k, rows and columns */
/* >                  k and -IPIV(k) were interchanged */
/* >                  in the matrix A(1:N,1:N). */
/* >                  If -IPIV(k) = k, no interchange occurred. */
/* >               2) If -IPIV(k-1) != k-1, rows and columns */
/* >                  k-1 and -IPIV(k-1) were interchanged */
/* >                  in the matrix A(1:N,1:N). */
/* >                  If -IPIV(k-1) = k-1, no interchange occurred. */
/* > */
/* >            c) In both cases a) and b), always ABS( IPIV(k) ) <= k. */
/* > */
/* >            d) NOTE: Any entry IPIV(k) is always NONZERO on output. */
/* > */
/* >          If UPLO = 'L', */
/* >          ( in factorization order, k increases from 1 to N ): */
/* >            a) A single positive entry IPIV(k) > 0 means: */
/* >               D(k,k) is a 1-by-1 diagonal block. */
/* >               If IPIV(k) != k, rows and columns k and IPIV(k) were */
/* >               interchanged in the matrix A(1:N,1:N). */
/* >               If IPIV(k) = k, no interchange occurred. */
/* > */
/* >            b) A pair of consecutive negative entries */
/* >               IPIV(k) < 0 and IPIV(k+1) < 0 means: */
/* >               D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* >               (NOTE: negative entries in IPIV appear ONLY in pairs). */
/* >               1) If -IPIV(k) != k, rows and columns */
/* >                  k and -IPIV(k) were interchanged */
/* >                  in the matrix A(1:N,1:N). */
/* >                  If -IPIV(k) = k, no interchange occurred. */
/* >               2) If -IPIV(k+1) != k+1, rows and columns */
/* >                  k-1 and -IPIV(k-1) were interchanged */
/* >                  in the matrix A(1:N,1:N). */
/* >                  If -IPIV(k+1) = k+1, no interchange occurred. */
/* > */
/* >            c) In both cases a) and b), always ABS( IPIV(k) ) >= k. */
/* > */
/* >            d) NOTE: Any entry IPIV(k) is always NONZERO on output. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* > */
/* >          < 0: If INFO = -k, the k-th argument had an illegal value */
/* > */
/* >          > 0: If INFO = k, the matrix A is singular, because: */
/* >                 If UPLO = 'U': column k in the upper */
/* >                 triangular part of A contains all zeros. */
/* >                 If UPLO = 'L': column k in the lower */
/* >                 triangular part of A contains all zeros. */
/* > */
/* >               Therefore D(k,k) is exactly zero, and superdiagonal */
/* >               elements of column k of U (or subdiagonal elements of */
/* >               column k of L ) are all zeros. The factorization has */
/* >               been completed, but the block diagonal matrix D is */
/* >               exactly singular, and division by zero will occur if */
/* >               it is used to solve a system of equations. */
/* > */
/* >               NOTE: INFO only stores the first occurrence of */
/* >               a singularity, any subsequent occurrence of singularity */
/* >               is not stored in INFO even though the factorization */
/* >               always completes. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleSYcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > TODO: put further details */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  December 2016,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* >  01-01-96 - Based on modifications by */
/* >    J. Lewis, Boeing Computer Services Company */
/* >    A. Petitet, Computer Science Dept., */
/* >                Univ. of Tenn., Knoxville abd , USA */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int dsytf2_rk__(char *uplo, integer *n, doublereal *a, 
	integer *lda, doublereal *e, integer *ipiv, integer *info, ftnlen 
	uplo_len)
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


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

/*     Test the input parameters. */

#line 289 "dsytf2_rk.f"
    /* Parameter adjustments */
#line 289 "dsytf2_rk.f"
    a_dim1 = *lda;
#line 289 "dsytf2_rk.f"
    a_offset = 1 + a_dim1;
#line 289 "dsytf2_rk.f"
    a -= a_offset;
#line 289 "dsytf2_rk.f"
    --e;
#line 289 "dsytf2_rk.f"
    --ipiv;
#line 289 "dsytf2_rk.f"

#line 289 "dsytf2_rk.f"
    /* Function Body */
#line 289 "dsytf2_rk.f"
    *info = 0;
#line 290 "dsytf2_rk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 291 "dsytf2_rk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 292 "dsytf2_rk.f"
	*info = -1;
#line 293 "dsytf2_rk.f"
    } else if (*n < 0) {
#line 294 "dsytf2_rk.f"
	*info = -2;
#line 295 "dsytf2_rk.f"
    } else if (*lda < max(1,*n)) {
#line 296 "dsytf2_rk.f"
	*info = -4;
#line 297 "dsytf2_rk.f"
    }
#line 298 "dsytf2_rk.f"
    if (*info != 0) {
#line 299 "dsytf2_rk.f"
	i__1 = -(*info);
#line 299 "dsytf2_rk.f"
	xerbla_("DSYTF2_RK", &i__1, (ftnlen)9);
#line 300 "dsytf2_rk.f"
	return 0;
#line 301 "dsytf2_rk.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 305 "dsytf2_rk.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 309 "dsytf2_rk.f"
    sfmin = dlamch_("S", (ftnlen)1);

#line 311 "dsytf2_rk.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        Initilize the first entry of array E, where superdiagonal */
/*        elements of D are stored */

#line 318 "dsytf2_rk.f"
	e[1] = 0.;

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 323 "dsytf2_rk.f"
	k = *n;
#line 324 "dsytf2_rk.f"
L10:

/*        If K < 1, exit from loop */

#line 328 "dsytf2_rk.f"
	if (k < 1) {
#line 328 "dsytf2_rk.f"
	    goto L34;
#line 328 "dsytf2_rk.f"
	}
#line 330 "dsytf2_rk.f"
	kstep = 1;
#line 331 "dsytf2_rk.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 336 "dsytf2_rk.f"
	absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 342 "dsytf2_rk.f"
	if (k > 1) {
#line 343 "dsytf2_rk.f"
	    i__1 = k - 1;
#line 343 "dsytf2_rk.f"
	    imax = idamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 344 "dsytf2_rk.f"
	    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
#line 345 "dsytf2_rk.f"
	} else {
#line 346 "dsytf2_rk.f"
	    colmax = 0.;
#line 347 "dsytf2_rk.f"
	}

#line 349 "dsytf2_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 353 "dsytf2_rk.f"
	    if (*info == 0) {
#line 353 "dsytf2_rk.f"
		*info = k;
#line 353 "dsytf2_rk.f"
	    }
#line 355 "dsytf2_rk.f"
	    kp = k;

/*           Set E( K ) to zero */

#line 359 "dsytf2_rk.f"
	    if (k > 1) {
#line 359 "dsytf2_rk.f"
		e[k] = 0.;
#line 359 "dsytf2_rk.f"
	    }

#line 362 "dsytf2_rk.f"
	} else {

/*           Test for interchange */

/*           Equivalent to testing for (used to handle NaN and Inf) */
/*           ABSAKK.GE.ALPHA*COLMAX */

#line 369 "dsytf2_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, */
/*              use 1-by-1 pivot block */

#line 374 "dsytf2_rk.f"
		kp = k;
#line 375 "dsytf2_rk.f"
	    } else {

#line 377 "dsytf2_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 381 "dsytf2_rk.f"
L12:

/*                 Begin pivot search loop body */

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 389 "dsytf2_rk.f"
		if (imax != k) {
#line 390 "dsytf2_rk.f"
		    i__1 = k - imax;
#line 390 "dsytf2_rk.f"
		    jmax = imax + idamax_(&i__1, &a[imax + (imax + 1) * 
			    a_dim1], lda);
#line 392 "dsytf2_rk.f"
		    rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
#line 393 "dsytf2_rk.f"
		} else {
#line 394 "dsytf2_rk.f"
		    rowmax = 0.;
#line 395 "dsytf2_rk.f"
		}

#line 397 "dsytf2_rk.f"
		if (imax > 1) {
#line 398 "dsytf2_rk.f"
		    i__1 = imax - 1;
#line 398 "dsytf2_rk.f"
		    itemp = idamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
#line 399 "dsytf2_rk.f"
		    dtemp = (d__1 = a[itemp + imax * a_dim1], abs(d__1));
#line 400 "dsytf2_rk.f"
		    if (dtemp > rowmax) {
#line 401 "dsytf2_rk.f"
			rowmax = dtemp;
#line 402 "dsytf2_rk.f"
			jmax = itemp;
#line 403 "dsytf2_rk.f"
		    }
#line 404 "dsytf2_rk.f"
		}

/*                 Equivalent to testing for (used to handle NaN and Inf) */
/*                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */

#line 409 "dsytf2_rk.f"
		if (! ((d__1 = a[imax + imax * a_dim1], abs(d__1)) < alpha * 
			rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 415 "dsytf2_rk.f"
		    kp = imax;
#line 416 "dsytf2_rk.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX .EQ. COLMAX, */
/*                 used to handle NaN and Inf */

#line 421 "dsytf2_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 426 "dsytf2_rk.f"
		    kp = imax;
#line 427 "dsytf2_rk.f"
		    kstep = 2;
#line 428 "dsytf2_rk.f"
		    done = TRUE_;
#line 429 "dsytf2_rk.f"
		} else {

/*                    Pivot NOT found, set variables and repeat */

#line 433 "dsytf2_rk.f"
		    p = imax;
#line 434 "dsytf2_rk.f"
		    colmax = rowmax;
#line 435 "dsytf2_rk.f"
		    imax = jmax;
#line 436 "dsytf2_rk.f"
		}

/*                 End pivot search loop body */

#line 440 "dsytf2_rk.f"
		if (! done) {
#line 440 "dsytf2_rk.f"
		    goto L12;
#line 440 "dsytf2_rk.f"
		}

#line 442 "dsytf2_rk.f"
	    }

/*           Swap TWO rows and TWO columns */

/*           First swap */

#line 448 "dsytf2_rk.f"
	    if (kstep == 2 && p != k) {

/*              Interchange rows and column K and P in the leading */
/*              submatrix A(1:k,1:k) if we have a 2-by-2 pivot */

#line 453 "dsytf2_rk.f"
		if (p > 1) {
#line 453 "dsytf2_rk.f"
		    i__1 = p - 1;
#line 453 "dsytf2_rk.f"
		    dswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 
			    1], &c__1);
#line 453 "dsytf2_rk.f"
		}
#line 455 "dsytf2_rk.f"
		if (p < k - 1) {
#line 455 "dsytf2_rk.f"
		    i__1 = k - p - 1;
#line 455 "dsytf2_rk.f"
		    dswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 
			    1) * a_dim1], lda);
#line 455 "dsytf2_rk.f"
		}
#line 458 "dsytf2_rk.f"
		t = a[k + k * a_dim1];
#line 459 "dsytf2_rk.f"
		a[k + k * a_dim1] = a[p + p * a_dim1];
#line 460 "dsytf2_rk.f"
		a[p + p * a_dim1] = t;

/*              Convert upper triangle of A into U form by applying */
/*              the interchanges in columns k+1:N. */

#line 465 "dsytf2_rk.f"
		if (k < *n) {
#line 465 "dsytf2_rk.f"
		    i__1 = *n - k;
#line 465 "dsytf2_rk.f"
		    dswap_(&i__1, &a[k + (k + 1) * a_dim1], lda, &a[p + (k + 
			    1) * a_dim1], lda);
#line 465 "dsytf2_rk.f"
		}

#line 468 "dsytf2_rk.f"
	    }

/*           Second swap */

#line 472 "dsytf2_rk.f"
	    kk = k - kstep + 1;
#line 473 "dsytf2_rk.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 478 "dsytf2_rk.f"
		if (kp > 1) {
#line 478 "dsytf2_rk.f"
		    i__1 = kp - 1;
#line 478 "dsytf2_rk.f"
		    dswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 478 "dsytf2_rk.f"
		}
#line 480 "dsytf2_rk.f"
		if (kk > 1 && kp < kk - 1) {
#line 480 "dsytf2_rk.f"
		    i__1 = kk - kp - 1;
#line 480 "dsytf2_rk.f"
		    dswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (
			    kp + 1) * a_dim1], lda);
#line 480 "dsytf2_rk.f"
		}
#line 483 "dsytf2_rk.f"
		t = a[kk + kk * a_dim1];
#line 484 "dsytf2_rk.f"
		a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
#line 485 "dsytf2_rk.f"
		a[kp + kp * a_dim1] = t;
#line 486 "dsytf2_rk.f"
		if (kstep == 2) {
#line 487 "dsytf2_rk.f"
		    t = a[k - 1 + k * a_dim1];
#line 488 "dsytf2_rk.f"
		    a[k - 1 + k * a_dim1] = a[kp + k * a_dim1];
#line 489 "dsytf2_rk.f"
		    a[kp + k * a_dim1] = t;
#line 490 "dsytf2_rk.f"
		}

/*              Convert upper triangle of A into U form by applying */
/*              the interchanges in columns k+1:N. */

#line 495 "dsytf2_rk.f"
		if (k < *n) {
#line 495 "dsytf2_rk.f"
		    i__1 = *n - k;
#line 495 "dsytf2_rk.f"
		    dswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 495 "dsytf2_rk.f"
		}

#line 499 "dsytf2_rk.f"
	    }

/*           Update the leading submatrix */

#line 503 "dsytf2_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

#line 511 "dsytf2_rk.f"
		if (k > 1) {

/*                 Perform a rank-1 update of A(1:k-1,1:k-1) and */
/*                 store U(k) in column k */

#line 516 "dsytf2_rk.f"
		    if ((d__1 = a[k + k * a_dim1], abs(d__1)) >= sfmin) {

/*                    Perform a rank-1 update of A(1:k-1,1:k-1) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*1/D(k)*W(k)**T */

#line 522 "dsytf2_rk.f"
			d11 = 1. / a[k + k * a_dim1];
#line 523 "dsytf2_rk.f"
			i__1 = k - 1;
#line 523 "dsytf2_rk.f"
			d__1 = -d11;
#line 523 "dsytf2_rk.f"
			dsyr_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);

/*                    Store U(k) in column k */

#line 527 "dsytf2_rk.f"
			i__1 = k - 1;
#line 527 "dsytf2_rk.f"
			dscal_(&i__1, &d11, &a[k * a_dim1 + 1], &c__1);
#line 528 "dsytf2_rk.f"
		    } else {

/*                    Store L(k) in column K */

#line 532 "dsytf2_rk.f"
			d11 = a[k + k * a_dim1];
#line 533 "dsytf2_rk.f"
			i__1 = k - 1;
#line 533 "dsytf2_rk.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 534 "dsytf2_rk.f"
			    a[ii + k * a_dim1] /= d11;
#line 535 "dsytf2_rk.f"
/* L16: */
#line 535 "dsytf2_rk.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 542 "dsytf2_rk.f"
			i__1 = k - 1;
#line 542 "dsytf2_rk.f"
			d__1 = -d11;
#line 542 "dsytf2_rk.f"
			dsyr_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);
#line 543 "dsytf2_rk.f"
		    }

/*                 Store the superdiagonal element of D in array E */

#line 547 "dsytf2_rk.f"
		    e[k] = 0.;

#line 549 "dsytf2_rk.f"
		}

#line 551 "dsytf2_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 567 "dsytf2_rk.f"
		if (k > 2) {

#line 569 "dsytf2_rk.f"
		    d12 = a[k - 1 + k * a_dim1];
#line 570 "dsytf2_rk.f"
		    d22 = a[k - 1 + (k - 1) * a_dim1] / d12;
#line 571 "dsytf2_rk.f"
		    d11 = a[k + k * a_dim1] / d12;
#line 572 "dsytf2_rk.f"
		    t = 1. / (d11 * d22 - 1.);

#line 574 "dsytf2_rk.f"
		    for (j = k - 2; j >= 1; --j) {

#line 576 "dsytf2_rk.f"
			wkm1 = t * (d11 * a[j + (k - 1) * a_dim1] - a[j + k * 
				a_dim1]);
#line 577 "dsytf2_rk.f"
			wk = t * (d22 * a[j + k * a_dim1] - a[j + (k - 1) * 
				a_dim1]);

#line 579 "dsytf2_rk.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 580 "dsytf2_rk.f"
			    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ 
				    + k * a_dim1] / d12 * wk - a[i__ + (k - 1)
				     * a_dim1] / d12 * wkm1;
#line 582 "dsytf2_rk.f"
/* L20: */
#line 582 "dsytf2_rk.f"
			}

/*                    Store U(k) and U(k-1) in cols k and k-1 for row J */

#line 586 "dsytf2_rk.f"
			a[j + k * a_dim1] = wk / d12;
#line 587 "dsytf2_rk.f"
			a[j + (k - 1) * a_dim1] = wkm1 / d12;

#line 589 "dsytf2_rk.f"
/* L30: */
#line 589 "dsytf2_rk.f"
		    }

#line 591 "dsytf2_rk.f"
		}

/*              Copy superdiagonal elements of D(K) to E(K) and */
/*              ZERO out superdiagonal entry of A */

#line 596 "dsytf2_rk.f"
		e[k] = a[k - 1 + k * a_dim1];
#line 597 "dsytf2_rk.f"
		e[k - 1] = 0.;
#line 598 "dsytf2_rk.f"
		a[k - 1 + k * a_dim1] = 0.;

#line 600 "dsytf2_rk.f"
	    }

/*           End column K is nonsingular */

#line 604 "dsytf2_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 608 "dsytf2_rk.f"
	if (kstep == 1) {
#line 609 "dsytf2_rk.f"
	    ipiv[k] = kp;
#line 610 "dsytf2_rk.f"
	} else {
#line 611 "dsytf2_rk.f"
	    ipiv[k] = -p;
#line 612 "dsytf2_rk.f"
	    ipiv[k - 1] = -kp;
#line 613 "dsytf2_rk.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 617 "dsytf2_rk.f"
	k -= kstep;
#line 618 "dsytf2_rk.f"
	goto L10;

#line 620 "dsytf2_rk.f"
L34:

#line 622 "dsytf2_rk.f"
	;
#line 622 "dsytf2_rk.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        Initilize the unused last entry of the subdiagonal array E. */

#line 628 "dsytf2_rk.f"
	e[*n] = 0.;

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 633 "dsytf2_rk.f"
	k = 1;
#line 634 "dsytf2_rk.f"
L40:

/*        If K > N, exit from loop */

#line 638 "dsytf2_rk.f"
	if (k > *n) {
#line 638 "dsytf2_rk.f"
	    goto L64;
#line 638 "dsytf2_rk.f"
	}
#line 640 "dsytf2_rk.f"
	kstep = 1;
#line 641 "dsytf2_rk.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 646 "dsytf2_rk.f"
	absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 652 "dsytf2_rk.f"
	if (k < *n) {
#line 653 "dsytf2_rk.f"
	    i__1 = *n - k;
#line 653 "dsytf2_rk.f"
	    imax = k + idamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 654 "dsytf2_rk.f"
	    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
#line 655 "dsytf2_rk.f"
	} else {
#line 656 "dsytf2_rk.f"
	    colmax = 0.;
#line 657 "dsytf2_rk.f"
	}

#line 659 "dsytf2_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 663 "dsytf2_rk.f"
	    if (*info == 0) {
#line 663 "dsytf2_rk.f"
		*info = k;
#line 663 "dsytf2_rk.f"
	    }
#line 665 "dsytf2_rk.f"
	    kp = k;

/*           Set E( K ) to zero */

#line 669 "dsytf2_rk.f"
	    if (k < *n) {
#line 669 "dsytf2_rk.f"
		e[k] = 0.;
#line 669 "dsytf2_rk.f"
	    }

#line 672 "dsytf2_rk.f"
	} else {

/*           Test for interchange */

/*           Equivalent to testing for (used to handle NaN and Inf) */
/*           ABSAKK.GE.ALPHA*COLMAX */

#line 679 "dsytf2_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 683 "dsytf2_rk.f"
		kp = k;

#line 685 "dsytf2_rk.f"
	    } else {

#line 687 "dsytf2_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 691 "dsytf2_rk.f"
L42:

/*                 Begin pivot search loop body */

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 699 "dsytf2_rk.f"
		if (imax != k) {
#line 700 "dsytf2_rk.f"
		    i__1 = imax - k;
#line 700 "dsytf2_rk.f"
		    jmax = k - 1 + idamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 701 "dsytf2_rk.f"
		    rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
#line 702 "dsytf2_rk.f"
		} else {
#line 703 "dsytf2_rk.f"
		    rowmax = 0.;
#line 704 "dsytf2_rk.f"
		}

#line 706 "dsytf2_rk.f"
		if (imax < *n) {
#line 707 "dsytf2_rk.f"
		    i__1 = *n - imax;
#line 707 "dsytf2_rk.f"
		    itemp = imax + idamax_(&i__1, &a[imax + 1 + imax * a_dim1]
			    , &c__1);
#line 709 "dsytf2_rk.f"
		    dtemp = (d__1 = a[itemp + imax * a_dim1], abs(d__1));
#line 710 "dsytf2_rk.f"
		    if (dtemp > rowmax) {
#line 711 "dsytf2_rk.f"
			rowmax = dtemp;
#line 712 "dsytf2_rk.f"
			jmax = itemp;
#line 713 "dsytf2_rk.f"
		    }
#line 714 "dsytf2_rk.f"
		}

/*                 Equivalent to testing for (used to handle NaN and Inf) */
/*                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */

#line 719 "dsytf2_rk.f"
		if (! ((d__1 = a[imax + imax * a_dim1], abs(d__1)) < alpha * 
			rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 725 "dsytf2_rk.f"
		    kp = imax;
#line 726 "dsytf2_rk.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX .EQ. COLMAX, */
/*                 used to handle NaN and Inf */

#line 731 "dsytf2_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 736 "dsytf2_rk.f"
		    kp = imax;
#line 737 "dsytf2_rk.f"
		    kstep = 2;
#line 738 "dsytf2_rk.f"
		    done = TRUE_;
#line 739 "dsytf2_rk.f"
		} else {

/*                    Pivot NOT found, set variables and repeat */

#line 743 "dsytf2_rk.f"
		    p = imax;
#line 744 "dsytf2_rk.f"
		    colmax = rowmax;
#line 745 "dsytf2_rk.f"
		    imax = jmax;
#line 746 "dsytf2_rk.f"
		}

/*                 End pivot search loop body */

#line 750 "dsytf2_rk.f"
		if (! done) {
#line 750 "dsytf2_rk.f"
		    goto L42;
#line 750 "dsytf2_rk.f"
		}

#line 752 "dsytf2_rk.f"
	    }

/*           Swap TWO rows and TWO columns */

/*           First swap */

#line 758 "dsytf2_rk.f"
	    if (kstep == 2 && p != k) {

/*              Interchange rows and column K and P in the trailing */
/*              submatrix A(k:n,k:n) if we have a 2-by-2 pivot */

#line 763 "dsytf2_rk.f"
		if (p < *n) {
#line 763 "dsytf2_rk.f"
		    i__1 = *n - p;
#line 763 "dsytf2_rk.f"
		    dswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p 
			    * a_dim1], &c__1);
#line 763 "dsytf2_rk.f"
		}
#line 765 "dsytf2_rk.f"
		if (p > k + 1) {
#line 765 "dsytf2_rk.f"
		    i__1 = p - k - 1;
#line 765 "dsytf2_rk.f"
		    dswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[p + (k + 
			    1) * a_dim1], lda);
#line 765 "dsytf2_rk.f"
		}
#line 767 "dsytf2_rk.f"
		t = a[k + k * a_dim1];
#line 768 "dsytf2_rk.f"
		a[k + k * a_dim1] = a[p + p * a_dim1];
#line 769 "dsytf2_rk.f"
		a[p + p * a_dim1] = t;

/*              Convert lower triangle of A into L form by applying */
/*              the interchanges in columns 1:k-1. */

#line 774 "dsytf2_rk.f"
		if (k > 1) {
#line 774 "dsytf2_rk.f"
		    i__1 = k - 1;
#line 774 "dsytf2_rk.f"
		    dswap_(&i__1, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 774 "dsytf2_rk.f"
		}

#line 777 "dsytf2_rk.f"
	    }

/*           Second swap */

#line 781 "dsytf2_rk.f"
	    kk = k + kstep - 1;
#line 782 "dsytf2_rk.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 787 "dsytf2_rk.f"
		if (kp < *n) {
#line 787 "dsytf2_rk.f"
		    i__1 = *n - kp;
#line 787 "dsytf2_rk.f"
		    dswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 787 "dsytf2_rk.f"
		}
#line 789 "dsytf2_rk.f"
		if (kk < *n && kp > kk + 1) {
#line 789 "dsytf2_rk.f"
		    i__1 = kp - kk - 1;
#line 789 "dsytf2_rk.f"
		    dswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (
			    kk + 1) * a_dim1], lda);
#line 789 "dsytf2_rk.f"
		}
#line 792 "dsytf2_rk.f"
		t = a[kk + kk * a_dim1];
#line 793 "dsytf2_rk.f"
		a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
#line 794 "dsytf2_rk.f"
		a[kp + kp * a_dim1] = t;
#line 795 "dsytf2_rk.f"
		if (kstep == 2) {
#line 796 "dsytf2_rk.f"
		    t = a[k + 1 + k * a_dim1];
#line 797 "dsytf2_rk.f"
		    a[k + 1 + k * a_dim1] = a[kp + k * a_dim1];
#line 798 "dsytf2_rk.f"
		    a[kp + k * a_dim1] = t;
#line 799 "dsytf2_rk.f"
		}

/*              Convert lower triangle of A into L form by applying */
/*              the interchanges in columns 1:k-1. */

#line 804 "dsytf2_rk.f"
		if (k > 1) {
#line 804 "dsytf2_rk.f"
		    i__1 = k - 1;
#line 804 "dsytf2_rk.f"
		    dswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 804 "dsytf2_rk.f"
		}

#line 807 "dsytf2_rk.f"
	    }

/*           Update the trailing submatrix */

#line 811 "dsytf2_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 819 "dsytf2_rk.f"
		if (k < *n) {

/*              Perform a rank-1 update of A(k+1:n,k+1:n) and */
/*              store L(k) in column k */

#line 824 "dsytf2_rk.f"
		    if ((d__1 = a[k + k * a_dim1], abs(d__1)) >= sfmin) {

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */

#line 830 "dsytf2_rk.f"
			d11 = 1. / a[k + k * a_dim1];
#line 831 "dsytf2_rk.f"
			i__1 = *n - k;
#line 831 "dsytf2_rk.f"
			d__1 = -d11;
#line 831 "dsytf2_rk.f"
			dsyr_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);

/*                    Store L(k) in column k */

#line 836 "dsytf2_rk.f"
			i__1 = *n - k;
#line 836 "dsytf2_rk.f"
			dscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
#line 837 "dsytf2_rk.f"
		    } else {

/*                    Store L(k) in column k */

#line 841 "dsytf2_rk.f"
			d11 = a[k + k * a_dim1];
#line 842 "dsytf2_rk.f"
			i__1 = *n;
#line 842 "dsytf2_rk.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 843 "dsytf2_rk.f"
			    a[ii + k * a_dim1] /= d11;
#line 844 "dsytf2_rk.f"
/* L46: */
#line 844 "dsytf2_rk.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 851 "dsytf2_rk.f"
			i__1 = *n - k;
#line 851 "dsytf2_rk.f"
			d__1 = -d11;
#line 851 "dsytf2_rk.f"
			dsyr_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);
#line 853 "dsytf2_rk.f"
		    }

/*                 Store the subdiagonal element of D in array E */

#line 857 "dsytf2_rk.f"
		    e[k] = 0.;

#line 859 "dsytf2_rk.f"
		}

#line 861 "dsytf2_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */


/*              Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T */
/*                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 878 "dsytf2_rk.f"
		if (k < *n - 1) {

#line 880 "dsytf2_rk.f"
		    d21 = a[k + 1 + k * a_dim1];
#line 881 "dsytf2_rk.f"
		    d11 = a[k + 1 + (k + 1) * a_dim1] / d21;
#line 882 "dsytf2_rk.f"
		    d22 = a[k + k * a_dim1] / d21;
#line 883 "dsytf2_rk.f"
		    t = 1. / (d11 * d22 - 1.);

#line 885 "dsytf2_rk.f"
		    i__1 = *n;
#line 885 "dsytf2_rk.f"
		    for (j = k + 2; j <= i__1; ++j) {

/*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */

#line 889 "dsytf2_rk.f"
			wk = t * (d11 * a[j + k * a_dim1] - a[j + (k + 1) * 
				a_dim1]);
#line 890 "dsytf2_rk.f"
			wkp1 = t * (d22 * a[j + (k + 1) * a_dim1] - a[j + k * 
				a_dim1]);

/*                    Perform a rank-2 update of A(k+2:n,k+2:n) */

#line 894 "dsytf2_rk.f"
			i__2 = *n;
#line 894 "dsytf2_rk.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 895 "dsytf2_rk.f"
			    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ 
				    + k * a_dim1] / d21 * wk - a[i__ + (k + 1)
				     * a_dim1] / d21 * wkp1;
#line 897 "dsytf2_rk.f"
/* L50: */
#line 897 "dsytf2_rk.f"
			}

/*                    Store L(k) and L(k+1) in cols k and k+1 for row J */

#line 901 "dsytf2_rk.f"
			a[j + k * a_dim1] = wk / d21;
#line 902 "dsytf2_rk.f"
			a[j + (k + 1) * a_dim1] = wkp1 / d21;

#line 904 "dsytf2_rk.f"
/* L60: */
#line 904 "dsytf2_rk.f"
		    }

#line 906 "dsytf2_rk.f"
		}

/*              Copy subdiagonal elements of D(K) to E(K) and */
/*              ZERO out subdiagonal entry of A */

#line 911 "dsytf2_rk.f"
		e[k] = a[k + 1 + k * a_dim1];
#line 912 "dsytf2_rk.f"
		e[k + 1] = 0.;
#line 913 "dsytf2_rk.f"
		a[k + 1 + k * a_dim1] = 0.;

#line 915 "dsytf2_rk.f"
	    }

/*           End column K is nonsingular */

#line 919 "dsytf2_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 923 "dsytf2_rk.f"
	if (kstep == 1) {
#line 924 "dsytf2_rk.f"
	    ipiv[k] = kp;
#line 925 "dsytf2_rk.f"
	} else {
#line 926 "dsytf2_rk.f"
	    ipiv[k] = -p;
#line 927 "dsytf2_rk.f"
	    ipiv[k + 1] = -kp;
#line 928 "dsytf2_rk.f"
	}

/*        Increase K and return to the start of the main loop */

#line 932 "dsytf2_rk.f"
	k += kstep;
#line 933 "dsytf2_rk.f"
	goto L40;

#line 935 "dsytf2_rk.f"
L64:

#line 937 "dsytf2_rk.f"
	;
#line 937 "dsytf2_rk.f"
    }

#line 939 "dsytf2_rk.f"
    return 0;

/*     End of DSYTF2_RK */

} /* dsytf2_rk__ */

