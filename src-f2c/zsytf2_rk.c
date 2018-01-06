#line 1 "zsytf2_rk.f"
/* zsytf2_rk.f -- translated by f2c (version 20100827).
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

#line 1 "zsytf2_rk.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZSYTF2_RK computes the factorization of a complex symmetric indefinite matrix using the bounded
 Bunch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTF2_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytf2_
rk.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytf2_
rk.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytf2_
rk.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), E ( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > ZSYTF2_RK computes the factorization of a complex symmetric matrix A */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          E is COMPLEX*16 array, dimension (N) */
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

/* > \ingroup complex16SYcomputational */

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
/* Subroutine */ int zsytf2_rk__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *e, integer *ipiv, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, p;
    static doublecomplex t, d11, d12, d21, d22;
    static integer ii, kk, kp;
    static doublecomplex wk, wkm1, wkp1;
    static logical done;
    static integer imax, jmax;
    extern /* Subroutine */ int zsyr_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dtemp, sfmin;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer itemp, kstep;
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal absakk;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal colmax;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static doublereal rowmax;


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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 298 "zsytf2_rk.f"
    /* Parameter adjustments */
#line 298 "zsytf2_rk.f"
    a_dim1 = *lda;
#line 298 "zsytf2_rk.f"
    a_offset = 1 + a_dim1;
#line 298 "zsytf2_rk.f"
    a -= a_offset;
#line 298 "zsytf2_rk.f"
    --e;
#line 298 "zsytf2_rk.f"
    --ipiv;
#line 298 "zsytf2_rk.f"

#line 298 "zsytf2_rk.f"
    /* Function Body */
#line 298 "zsytf2_rk.f"
    *info = 0;
#line 299 "zsytf2_rk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 300 "zsytf2_rk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 301 "zsytf2_rk.f"
	*info = -1;
#line 302 "zsytf2_rk.f"
    } else if (*n < 0) {
#line 303 "zsytf2_rk.f"
	*info = -2;
#line 304 "zsytf2_rk.f"
    } else if (*lda < max(1,*n)) {
#line 305 "zsytf2_rk.f"
	*info = -4;
#line 306 "zsytf2_rk.f"
    }
#line 307 "zsytf2_rk.f"
    if (*info != 0) {
#line 308 "zsytf2_rk.f"
	i__1 = -(*info);
#line 308 "zsytf2_rk.f"
	xerbla_("ZSYTF2_RK", &i__1, (ftnlen)9);
#line 309 "zsytf2_rk.f"
	return 0;
#line 310 "zsytf2_rk.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 314 "zsytf2_rk.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 318 "zsytf2_rk.f"
    sfmin = dlamch_("S", (ftnlen)1);

#line 320 "zsytf2_rk.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        Initilize the first entry of array E, where superdiagonal */
/*        elements of D are stored */

#line 327 "zsytf2_rk.f"
	e[1].r = 0., e[1].i = 0.;

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 332 "zsytf2_rk.f"
	k = *n;
#line 333 "zsytf2_rk.f"
L10:

/*        If K < 1, exit from loop */

#line 337 "zsytf2_rk.f"
	if (k < 1) {
#line 337 "zsytf2_rk.f"
	    goto L34;
#line 337 "zsytf2_rk.f"
	}
#line 339 "zsytf2_rk.f"
	kstep = 1;
#line 340 "zsytf2_rk.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 345 "zsytf2_rk.f"
	i__1 = k + k * a_dim1;
#line 345 "zsytf2_rk.f"
	absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * 
		a_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 351 "zsytf2_rk.f"
	if (k > 1) {
#line 352 "zsytf2_rk.f"
	    i__1 = k - 1;
#line 352 "zsytf2_rk.f"
	    imax = izamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 353 "zsytf2_rk.f"
	    i__1 = imax + k * a_dim1;
#line 353 "zsytf2_rk.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 354 "zsytf2_rk.f"
	} else {
#line 355 "zsytf2_rk.f"
	    colmax = 0.;
#line 356 "zsytf2_rk.f"
	}

#line 358 "zsytf2_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 362 "zsytf2_rk.f"
	    if (*info == 0) {
#line 362 "zsytf2_rk.f"
		*info = k;
#line 362 "zsytf2_rk.f"
	    }
#line 364 "zsytf2_rk.f"
	    kp = k;

/*           Set E( K ) to zero */

#line 368 "zsytf2_rk.f"
	    if (k > 1) {
#line 368 "zsytf2_rk.f"
		i__1 = k;
#line 368 "zsytf2_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 368 "zsytf2_rk.f"
	    }

#line 371 "zsytf2_rk.f"
	} else {

/*           Test for interchange */

/*           Equivalent to testing for (used to handle NaN and Inf) */
/*           ABSAKK.GE.ALPHA*COLMAX */

#line 378 "zsytf2_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, */
/*              use 1-by-1 pivot block */

#line 383 "zsytf2_rk.f"
		kp = k;
#line 384 "zsytf2_rk.f"
	    } else {

#line 386 "zsytf2_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 390 "zsytf2_rk.f"
L12:

/*                 Begin pivot search loop body */

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 398 "zsytf2_rk.f"
		if (imax != k) {
#line 399 "zsytf2_rk.f"
		    i__1 = k - imax;
#line 399 "zsytf2_rk.f"
		    jmax = imax + izamax_(&i__1, &a[imax + (imax + 1) * 
			    a_dim1], lda);
#line 401 "zsytf2_rk.f"
		    i__1 = imax + jmax * a_dim1;
#line 401 "zsytf2_rk.f"
		    rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    a[imax + jmax * a_dim1]), abs(d__2));
#line 402 "zsytf2_rk.f"
		} else {
#line 403 "zsytf2_rk.f"
		    rowmax = 0.;
#line 404 "zsytf2_rk.f"
		}

#line 406 "zsytf2_rk.f"
		if (imax > 1) {
#line 407 "zsytf2_rk.f"
		    i__1 = imax - 1;
#line 407 "zsytf2_rk.f"
		    itemp = izamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
#line 408 "zsytf2_rk.f"
		    i__1 = itemp + imax * a_dim1;
#line 408 "zsytf2_rk.f"
		    dtemp = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    itemp + imax * a_dim1]), abs(d__2));
#line 409 "zsytf2_rk.f"
		    if (dtemp > rowmax) {
#line 410 "zsytf2_rk.f"
			rowmax = dtemp;
#line 411 "zsytf2_rk.f"
			jmax = itemp;
#line 412 "zsytf2_rk.f"
		    }
#line 413 "zsytf2_rk.f"
		}

/*                 Equivalent to testing for (used to handle NaN and Inf) */
/*                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */

#line 418 "zsytf2_rk.f"
		i__1 = imax + imax * a_dim1;
#line 418 "zsytf2_rk.f"
		if (! ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax 
			+ imax * a_dim1]), abs(d__2)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 424 "zsytf2_rk.f"
		    kp = imax;
#line 425 "zsytf2_rk.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX .EQ. COLMAX, */
/*                 used to handle NaN and Inf */

#line 430 "zsytf2_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 435 "zsytf2_rk.f"
		    kp = imax;
#line 436 "zsytf2_rk.f"
		    kstep = 2;
#line 437 "zsytf2_rk.f"
		    done = TRUE_;
#line 438 "zsytf2_rk.f"
		} else {

/*                    Pivot NOT found, set variables and repeat */

#line 442 "zsytf2_rk.f"
		    p = imax;
#line 443 "zsytf2_rk.f"
		    colmax = rowmax;
#line 444 "zsytf2_rk.f"
		    imax = jmax;
#line 445 "zsytf2_rk.f"
		}

/*                 End pivot search loop body */

#line 449 "zsytf2_rk.f"
		if (! done) {
#line 449 "zsytf2_rk.f"
		    goto L12;
#line 449 "zsytf2_rk.f"
		}

#line 451 "zsytf2_rk.f"
	    }

/*           Swap TWO rows and TWO columns */

/*           First swap */

#line 457 "zsytf2_rk.f"
	    if (kstep == 2 && p != k) {

/*              Interchange rows and column K and P in the leading */
/*              submatrix A(1:k,1:k) if we have a 2-by-2 pivot */

#line 462 "zsytf2_rk.f"
		if (p > 1) {
#line 462 "zsytf2_rk.f"
		    i__1 = p - 1;
#line 462 "zsytf2_rk.f"
		    zswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 
			    1], &c__1);
#line 462 "zsytf2_rk.f"
		}
#line 464 "zsytf2_rk.f"
		if (p < k - 1) {
#line 464 "zsytf2_rk.f"
		    i__1 = k - p - 1;
#line 464 "zsytf2_rk.f"
		    zswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 
			    1) * a_dim1], lda);
#line 464 "zsytf2_rk.f"
		}
#line 467 "zsytf2_rk.f"
		i__1 = k + k * a_dim1;
#line 467 "zsytf2_rk.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 468 "zsytf2_rk.f"
		i__1 = k + k * a_dim1;
#line 468 "zsytf2_rk.f"
		i__2 = p + p * a_dim1;
#line 468 "zsytf2_rk.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 469 "zsytf2_rk.f"
		i__1 = p + p * a_dim1;
#line 469 "zsytf2_rk.f"
		a[i__1].r = t.r, a[i__1].i = t.i;

/*              Convert upper triangle of A into U form by applying */
/*              the interchanges in columns k+1:N. */

#line 474 "zsytf2_rk.f"
		if (k < *n) {
#line 474 "zsytf2_rk.f"
		    i__1 = *n - k;
#line 474 "zsytf2_rk.f"
		    zswap_(&i__1, &a[k + (k + 1) * a_dim1], lda, &a[p + (k + 
			    1) * a_dim1], lda);
#line 474 "zsytf2_rk.f"
		}

#line 477 "zsytf2_rk.f"
	    }

/*           Second swap */

#line 481 "zsytf2_rk.f"
	    kk = k - kstep + 1;
#line 482 "zsytf2_rk.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 487 "zsytf2_rk.f"
		if (kp > 1) {
#line 487 "zsytf2_rk.f"
		    i__1 = kp - 1;
#line 487 "zsytf2_rk.f"
		    zswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 487 "zsytf2_rk.f"
		}
#line 489 "zsytf2_rk.f"
		if (kk > 1 && kp < kk - 1) {
#line 489 "zsytf2_rk.f"
		    i__1 = kk - kp - 1;
#line 489 "zsytf2_rk.f"
		    zswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (
			    kp + 1) * a_dim1], lda);
#line 489 "zsytf2_rk.f"
		}
#line 492 "zsytf2_rk.f"
		i__1 = kk + kk * a_dim1;
#line 492 "zsytf2_rk.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 493 "zsytf2_rk.f"
		i__1 = kk + kk * a_dim1;
#line 493 "zsytf2_rk.f"
		i__2 = kp + kp * a_dim1;
#line 493 "zsytf2_rk.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 494 "zsytf2_rk.f"
		i__1 = kp + kp * a_dim1;
#line 494 "zsytf2_rk.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 495 "zsytf2_rk.f"
		if (kstep == 2) {
#line 496 "zsytf2_rk.f"
		    i__1 = k - 1 + k * a_dim1;
#line 496 "zsytf2_rk.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 497 "zsytf2_rk.f"
		    i__1 = k - 1 + k * a_dim1;
#line 497 "zsytf2_rk.f"
		    i__2 = kp + k * a_dim1;
#line 497 "zsytf2_rk.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 498 "zsytf2_rk.f"
		    i__1 = kp + k * a_dim1;
#line 498 "zsytf2_rk.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 499 "zsytf2_rk.f"
		}

/*              Convert upper triangle of A into U form by applying */
/*              the interchanges in columns k+1:N. */

#line 504 "zsytf2_rk.f"
		if (k < *n) {
#line 504 "zsytf2_rk.f"
		    i__1 = *n - k;
#line 504 "zsytf2_rk.f"
		    zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 504 "zsytf2_rk.f"
		}

#line 508 "zsytf2_rk.f"
	    }

/*           Update the leading submatrix */

#line 512 "zsytf2_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

#line 520 "zsytf2_rk.f"
		if (k > 1) {

/*                 Perform a rank-1 update of A(1:k-1,1:k-1) and */
/*                 store U(k) in column k */

#line 525 "zsytf2_rk.f"
		    i__1 = k + k * a_dim1;
#line 525 "zsytf2_rk.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + 
			    k * a_dim1]), abs(d__2)) >= sfmin) {

/*                    Perform a rank-1 update of A(1:k-1,1:k-1) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*1/D(k)*W(k)**T */

#line 531 "zsytf2_rk.f"
			z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 531 "zsytf2_rk.f"
			d11.r = z__1.r, d11.i = z__1.i;
#line 532 "zsytf2_rk.f"
			i__1 = k - 1;
#line 532 "zsytf2_rk.f"
			z__1.r = -d11.r, z__1.i = -d11.i;
#line 532 "zsytf2_rk.f"
			zsyr_(uplo, &i__1, &z__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);

/*                    Store U(k) in column k */

#line 536 "zsytf2_rk.f"
			i__1 = k - 1;
#line 536 "zsytf2_rk.f"
			zscal_(&i__1, &d11, &a[k * a_dim1 + 1], &c__1);
#line 537 "zsytf2_rk.f"
		    } else {

/*                    Store L(k) in column K */

#line 541 "zsytf2_rk.f"
			i__1 = k + k * a_dim1;
#line 541 "zsytf2_rk.f"
			d11.r = a[i__1].r, d11.i = a[i__1].i;
#line 542 "zsytf2_rk.f"
			i__1 = k - 1;
#line 542 "zsytf2_rk.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 543 "zsytf2_rk.f"
			    i__2 = ii + k * a_dim1;
#line 543 "zsytf2_rk.f"
			    z_div(&z__1, &a[ii + k * a_dim1], &d11);
#line 543 "zsytf2_rk.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 544 "zsytf2_rk.f"
/* L16: */
#line 544 "zsytf2_rk.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 551 "zsytf2_rk.f"
			i__1 = k - 1;
#line 551 "zsytf2_rk.f"
			z__1.r = -d11.r, z__1.i = -d11.i;
#line 551 "zsytf2_rk.f"
			zsyr_(uplo, &i__1, &z__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);
#line 552 "zsytf2_rk.f"
		    }

/*                 Store the superdiagonal element of D in array E */

#line 556 "zsytf2_rk.f"
		    i__1 = k;
#line 556 "zsytf2_rk.f"
		    e[i__1].r = 0., e[i__1].i = 0.;

#line 558 "zsytf2_rk.f"
		}

#line 560 "zsytf2_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 576 "zsytf2_rk.f"
		if (k > 2) {

#line 578 "zsytf2_rk.f"
		    i__1 = k - 1 + k * a_dim1;
#line 578 "zsytf2_rk.f"
		    d12.r = a[i__1].r, d12.i = a[i__1].i;
#line 579 "zsytf2_rk.f"
		    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &d12);
#line 579 "zsytf2_rk.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 580 "zsytf2_rk.f"
		    z_div(&z__1, &a[k + k * a_dim1], &d12);
#line 580 "zsytf2_rk.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 581 "zsytf2_rk.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 581 "zsytf2_rk.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 581 "zsytf2_rk.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 581 "zsytf2_rk.f"
		    t.r = z__1.r, t.i = z__1.i;

#line 583 "zsytf2_rk.f"
		    for (j = k - 2; j >= 1; --j) {

#line 585 "zsytf2_rk.f"
			i__1 = j + (k - 1) * a_dim1;
#line 585 "zsytf2_rk.f"
			z__3.r = d11.r * a[i__1].r - d11.i * a[i__1].i, 
				z__3.i = d11.r * a[i__1].i + d11.i * a[i__1]
				.r;
#line 585 "zsytf2_rk.f"
			i__2 = j + k * a_dim1;
#line 585 "zsytf2_rk.f"
			z__2.r = z__3.r - a[i__2].r, z__2.i = z__3.i - a[i__2]
				.i;
#line 585 "zsytf2_rk.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 585 "zsytf2_rk.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 586 "zsytf2_rk.f"
			i__1 = j + k * a_dim1;
#line 586 "zsytf2_rk.f"
			z__3.r = d22.r * a[i__1].r - d22.i * a[i__1].i, 
				z__3.i = d22.r * a[i__1].i + d22.i * a[i__1]
				.r;
#line 586 "zsytf2_rk.f"
			i__2 = j + (k - 1) * a_dim1;
#line 586 "zsytf2_rk.f"
			z__2.r = z__3.r - a[i__2].r, z__2.i = z__3.i - a[i__2]
				.i;
#line 586 "zsytf2_rk.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 586 "zsytf2_rk.f"
			wk.r = z__1.r, wk.i = z__1.i;

#line 588 "zsytf2_rk.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 589 "zsytf2_rk.f"
			    i__1 = i__ + j * a_dim1;
#line 589 "zsytf2_rk.f"
			    i__2 = i__ + j * a_dim1;
#line 589 "zsytf2_rk.f"
			    z_div(&z__4, &a[i__ + k * a_dim1], &d12);
#line 589 "zsytf2_rk.f"
			    z__3.r = z__4.r * wk.r - z__4.i * wk.i, z__3.i = 
				    z__4.r * wk.i + z__4.i * wk.r;
#line 589 "zsytf2_rk.f"
			    z__2.r = a[i__2].r - z__3.r, z__2.i = a[i__2].i - 
				    z__3.i;
#line 589 "zsytf2_rk.f"
			    z_div(&z__6, &a[i__ + (k - 1) * a_dim1], &d12);
#line 589 "zsytf2_rk.f"
			    z__5.r = z__6.r * wkm1.r - z__6.i * wkm1.i, 
				    z__5.i = z__6.r * wkm1.i + z__6.i * 
				    wkm1.r;
#line 589 "zsytf2_rk.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 589 "zsytf2_rk.f"
			    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 591 "zsytf2_rk.f"
/* L20: */
#line 591 "zsytf2_rk.f"
			}

/*                    Store U(k) and U(k-1) in cols k and k-1 for row J */

#line 595 "zsytf2_rk.f"
			i__1 = j + k * a_dim1;
#line 595 "zsytf2_rk.f"
			z_div(&z__1, &wk, &d12);
#line 595 "zsytf2_rk.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 596 "zsytf2_rk.f"
			i__1 = j + (k - 1) * a_dim1;
#line 596 "zsytf2_rk.f"
			z_div(&z__1, &wkm1, &d12);
#line 596 "zsytf2_rk.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 598 "zsytf2_rk.f"
/* L30: */
#line 598 "zsytf2_rk.f"
		    }

#line 600 "zsytf2_rk.f"
		}

/*              Copy superdiagonal elements of D(K) to E(K) and */
/*              ZERO out superdiagonal entry of A */

#line 605 "zsytf2_rk.f"
		i__1 = k;
#line 605 "zsytf2_rk.f"
		i__2 = k - 1 + k * a_dim1;
#line 605 "zsytf2_rk.f"
		e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 606 "zsytf2_rk.f"
		i__1 = k - 1;
#line 606 "zsytf2_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 607 "zsytf2_rk.f"
		i__1 = k - 1 + k * a_dim1;
#line 607 "zsytf2_rk.f"
		a[i__1].r = 0., a[i__1].i = 0.;

#line 609 "zsytf2_rk.f"
	    }

/*           End column K is nonsingular */

#line 613 "zsytf2_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 617 "zsytf2_rk.f"
	if (kstep == 1) {
#line 618 "zsytf2_rk.f"
	    ipiv[k] = kp;
#line 619 "zsytf2_rk.f"
	} else {
#line 620 "zsytf2_rk.f"
	    ipiv[k] = -p;
#line 621 "zsytf2_rk.f"
	    ipiv[k - 1] = -kp;
#line 622 "zsytf2_rk.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 626 "zsytf2_rk.f"
	k -= kstep;
#line 627 "zsytf2_rk.f"
	goto L10;

#line 629 "zsytf2_rk.f"
L34:

#line 631 "zsytf2_rk.f"
	;
#line 631 "zsytf2_rk.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        Initilize the unused last entry of the subdiagonal array E. */

#line 637 "zsytf2_rk.f"
	i__1 = *n;
#line 637 "zsytf2_rk.f"
	e[i__1].r = 0., e[i__1].i = 0.;

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 642 "zsytf2_rk.f"
	k = 1;
#line 643 "zsytf2_rk.f"
L40:

/*        If K > N, exit from loop */

#line 647 "zsytf2_rk.f"
	if (k > *n) {
#line 647 "zsytf2_rk.f"
	    goto L64;
#line 647 "zsytf2_rk.f"
	}
#line 649 "zsytf2_rk.f"
	kstep = 1;
#line 650 "zsytf2_rk.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 655 "zsytf2_rk.f"
	i__1 = k + k * a_dim1;
#line 655 "zsytf2_rk.f"
	absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * 
		a_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 661 "zsytf2_rk.f"
	if (k < *n) {
#line 662 "zsytf2_rk.f"
	    i__1 = *n - k;
#line 662 "zsytf2_rk.f"
	    imax = k + izamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 663 "zsytf2_rk.f"
	    i__1 = imax + k * a_dim1;
#line 663 "zsytf2_rk.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 664 "zsytf2_rk.f"
	} else {
#line 665 "zsytf2_rk.f"
	    colmax = 0.;
#line 666 "zsytf2_rk.f"
	}

#line 668 "zsytf2_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 672 "zsytf2_rk.f"
	    if (*info == 0) {
#line 672 "zsytf2_rk.f"
		*info = k;
#line 672 "zsytf2_rk.f"
	    }
#line 674 "zsytf2_rk.f"
	    kp = k;

/*           Set E( K ) to zero */

#line 678 "zsytf2_rk.f"
	    if (k < *n) {
#line 678 "zsytf2_rk.f"
		i__1 = k;
#line 678 "zsytf2_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 678 "zsytf2_rk.f"
	    }

#line 681 "zsytf2_rk.f"
	} else {

/*           Test for interchange */

/*           Equivalent to testing for (used to handle NaN and Inf) */
/*           ABSAKK.GE.ALPHA*COLMAX */

#line 688 "zsytf2_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 692 "zsytf2_rk.f"
		kp = k;

#line 694 "zsytf2_rk.f"
	    } else {

#line 696 "zsytf2_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 700 "zsytf2_rk.f"
L42:

/*                 Begin pivot search loop body */

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 708 "zsytf2_rk.f"
		if (imax != k) {
#line 709 "zsytf2_rk.f"
		    i__1 = imax - k;
#line 709 "zsytf2_rk.f"
		    jmax = k - 1 + izamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 710 "zsytf2_rk.f"
		    i__1 = imax + jmax * a_dim1;
#line 710 "zsytf2_rk.f"
		    rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    a[imax + jmax * a_dim1]), abs(d__2));
#line 711 "zsytf2_rk.f"
		} else {
#line 712 "zsytf2_rk.f"
		    rowmax = 0.;
#line 713 "zsytf2_rk.f"
		}

#line 715 "zsytf2_rk.f"
		if (imax < *n) {
#line 716 "zsytf2_rk.f"
		    i__1 = *n - imax;
#line 716 "zsytf2_rk.f"
		    itemp = imax + izamax_(&i__1, &a[imax + 1 + imax * a_dim1]
			    , &c__1);
#line 718 "zsytf2_rk.f"
		    i__1 = itemp + imax * a_dim1;
#line 718 "zsytf2_rk.f"
		    dtemp = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    itemp + imax * a_dim1]), abs(d__2));
#line 719 "zsytf2_rk.f"
		    if (dtemp > rowmax) {
#line 720 "zsytf2_rk.f"
			rowmax = dtemp;
#line 721 "zsytf2_rk.f"
			jmax = itemp;
#line 722 "zsytf2_rk.f"
		    }
#line 723 "zsytf2_rk.f"
		}

/*                 Equivalent to testing for (used to handle NaN and Inf) */
/*                 ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */

#line 728 "zsytf2_rk.f"
		i__1 = imax + imax * a_dim1;
#line 728 "zsytf2_rk.f"
		if (! ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax 
			+ imax * a_dim1]), abs(d__2)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 734 "zsytf2_rk.f"
		    kp = imax;
#line 735 "zsytf2_rk.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX .EQ. COLMAX, */
/*                 used to handle NaN and Inf */

#line 740 "zsytf2_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 745 "zsytf2_rk.f"
		    kp = imax;
#line 746 "zsytf2_rk.f"
		    kstep = 2;
#line 747 "zsytf2_rk.f"
		    done = TRUE_;
#line 748 "zsytf2_rk.f"
		} else {

/*                    Pivot NOT found, set variables and repeat */

#line 752 "zsytf2_rk.f"
		    p = imax;
#line 753 "zsytf2_rk.f"
		    colmax = rowmax;
#line 754 "zsytf2_rk.f"
		    imax = jmax;
#line 755 "zsytf2_rk.f"
		}

/*                 End pivot search loop body */

#line 759 "zsytf2_rk.f"
		if (! done) {
#line 759 "zsytf2_rk.f"
		    goto L42;
#line 759 "zsytf2_rk.f"
		}

#line 761 "zsytf2_rk.f"
	    }

/*           Swap TWO rows and TWO columns */

/*           First swap */

#line 767 "zsytf2_rk.f"
	    if (kstep == 2 && p != k) {

/*              Interchange rows and column K and P in the trailing */
/*              submatrix A(k:n,k:n) if we have a 2-by-2 pivot */

#line 772 "zsytf2_rk.f"
		if (p < *n) {
#line 772 "zsytf2_rk.f"
		    i__1 = *n - p;
#line 772 "zsytf2_rk.f"
		    zswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p 
			    * a_dim1], &c__1);
#line 772 "zsytf2_rk.f"
		}
#line 774 "zsytf2_rk.f"
		if (p > k + 1) {
#line 774 "zsytf2_rk.f"
		    i__1 = p - k - 1;
#line 774 "zsytf2_rk.f"
		    zswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[p + (k + 
			    1) * a_dim1], lda);
#line 774 "zsytf2_rk.f"
		}
#line 776 "zsytf2_rk.f"
		i__1 = k + k * a_dim1;
#line 776 "zsytf2_rk.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 777 "zsytf2_rk.f"
		i__1 = k + k * a_dim1;
#line 777 "zsytf2_rk.f"
		i__2 = p + p * a_dim1;
#line 777 "zsytf2_rk.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 778 "zsytf2_rk.f"
		i__1 = p + p * a_dim1;
#line 778 "zsytf2_rk.f"
		a[i__1].r = t.r, a[i__1].i = t.i;

/*              Convert lower triangle of A into L form by applying */
/*              the interchanges in columns 1:k-1. */

#line 783 "zsytf2_rk.f"
		if (k > 1) {
#line 783 "zsytf2_rk.f"
		    i__1 = k - 1;
#line 783 "zsytf2_rk.f"
		    zswap_(&i__1, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 783 "zsytf2_rk.f"
		}

#line 786 "zsytf2_rk.f"
	    }

/*           Second swap */

#line 790 "zsytf2_rk.f"
	    kk = k + kstep - 1;
#line 791 "zsytf2_rk.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 796 "zsytf2_rk.f"
		if (kp < *n) {
#line 796 "zsytf2_rk.f"
		    i__1 = *n - kp;
#line 796 "zsytf2_rk.f"
		    zswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 796 "zsytf2_rk.f"
		}
#line 798 "zsytf2_rk.f"
		if (kk < *n && kp > kk + 1) {
#line 798 "zsytf2_rk.f"
		    i__1 = kp - kk - 1;
#line 798 "zsytf2_rk.f"
		    zswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (
			    kk + 1) * a_dim1], lda);
#line 798 "zsytf2_rk.f"
		}
#line 801 "zsytf2_rk.f"
		i__1 = kk + kk * a_dim1;
#line 801 "zsytf2_rk.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 802 "zsytf2_rk.f"
		i__1 = kk + kk * a_dim1;
#line 802 "zsytf2_rk.f"
		i__2 = kp + kp * a_dim1;
#line 802 "zsytf2_rk.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 803 "zsytf2_rk.f"
		i__1 = kp + kp * a_dim1;
#line 803 "zsytf2_rk.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 804 "zsytf2_rk.f"
		if (kstep == 2) {
#line 805 "zsytf2_rk.f"
		    i__1 = k + 1 + k * a_dim1;
#line 805 "zsytf2_rk.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 806 "zsytf2_rk.f"
		    i__1 = k + 1 + k * a_dim1;
#line 806 "zsytf2_rk.f"
		    i__2 = kp + k * a_dim1;
#line 806 "zsytf2_rk.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 807 "zsytf2_rk.f"
		    i__1 = kp + k * a_dim1;
#line 807 "zsytf2_rk.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 808 "zsytf2_rk.f"
		}

/*              Convert lower triangle of A into L form by applying */
/*              the interchanges in columns 1:k-1. */

#line 813 "zsytf2_rk.f"
		if (k > 1) {
#line 813 "zsytf2_rk.f"
		    i__1 = k - 1;
#line 813 "zsytf2_rk.f"
		    zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 813 "zsytf2_rk.f"
		}

#line 816 "zsytf2_rk.f"
	    }

/*           Update the trailing submatrix */

#line 820 "zsytf2_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 828 "zsytf2_rk.f"
		if (k < *n) {

/*              Perform a rank-1 update of A(k+1:n,k+1:n) and */
/*              store L(k) in column k */

#line 833 "zsytf2_rk.f"
		    i__1 = k + k * a_dim1;
#line 833 "zsytf2_rk.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + 
			    k * a_dim1]), abs(d__2)) >= sfmin) {

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */

#line 839 "zsytf2_rk.f"
			z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 839 "zsytf2_rk.f"
			d11.r = z__1.r, d11.i = z__1.i;
#line 840 "zsytf2_rk.f"
			i__1 = *n - k;
#line 840 "zsytf2_rk.f"
			z__1.r = -d11.r, z__1.i = -d11.i;
#line 840 "zsytf2_rk.f"
			zsyr_(uplo, &i__1, &z__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);

/*                    Store L(k) in column k */

#line 845 "zsytf2_rk.f"
			i__1 = *n - k;
#line 845 "zsytf2_rk.f"
			zscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
#line 846 "zsytf2_rk.f"
		    } else {

/*                    Store L(k) in column k */

#line 850 "zsytf2_rk.f"
			i__1 = k + k * a_dim1;
#line 850 "zsytf2_rk.f"
			d11.r = a[i__1].r, d11.i = a[i__1].i;
#line 851 "zsytf2_rk.f"
			i__1 = *n;
#line 851 "zsytf2_rk.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 852 "zsytf2_rk.f"
			    i__2 = ii + k * a_dim1;
#line 852 "zsytf2_rk.f"
			    z_div(&z__1, &a[ii + k * a_dim1], &d11);
#line 852 "zsytf2_rk.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 853 "zsytf2_rk.f"
/* L46: */
#line 853 "zsytf2_rk.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 860 "zsytf2_rk.f"
			i__1 = *n - k;
#line 860 "zsytf2_rk.f"
			z__1.r = -d11.r, z__1.i = -d11.i;
#line 860 "zsytf2_rk.f"
			zsyr_(uplo, &i__1, &z__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);
#line 862 "zsytf2_rk.f"
		    }

/*                 Store the subdiagonal element of D in array E */

#line 866 "zsytf2_rk.f"
		    i__1 = k;
#line 866 "zsytf2_rk.f"
		    e[i__1].r = 0., e[i__1].i = 0.;

#line 868 "zsytf2_rk.f"
		}

#line 870 "zsytf2_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */


/*              Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T */
/*                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 887 "zsytf2_rk.f"
		if (k < *n - 1) {

#line 889 "zsytf2_rk.f"
		    i__1 = k + 1 + k * a_dim1;
#line 889 "zsytf2_rk.f"
		    d21.r = a[i__1].r, d21.i = a[i__1].i;
#line 890 "zsytf2_rk.f"
		    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &d21);
#line 890 "zsytf2_rk.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 891 "zsytf2_rk.f"
		    z_div(&z__1, &a[k + k * a_dim1], &d21);
#line 891 "zsytf2_rk.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 892 "zsytf2_rk.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 892 "zsytf2_rk.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 892 "zsytf2_rk.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 892 "zsytf2_rk.f"
		    t.r = z__1.r, t.i = z__1.i;

#line 894 "zsytf2_rk.f"
		    i__1 = *n;
#line 894 "zsytf2_rk.f"
		    for (j = k + 2; j <= i__1; ++j) {

/*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */

#line 898 "zsytf2_rk.f"
			i__2 = j + k * a_dim1;
#line 898 "zsytf2_rk.f"
			z__3.r = d11.r * a[i__2].r - d11.i * a[i__2].i, 
				z__3.i = d11.r * a[i__2].i + d11.i * a[i__2]
				.r;
#line 898 "zsytf2_rk.f"
			i__3 = j + (k + 1) * a_dim1;
#line 898 "zsytf2_rk.f"
			z__2.r = z__3.r - a[i__3].r, z__2.i = z__3.i - a[i__3]
				.i;
#line 898 "zsytf2_rk.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 898 "zsytf2_rk.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 899 "zsytf2_rk.f"
			i__2 = j + (k + 1) * a_dim1;
#line 899 "zsytf2_rk.f"
			z__3.r = d22.r * a[i__2].r - d22.i * a[i__2].i, 
				z__3.i = d22.r * a[i__2].i + d22.i * a[i__2]
				.r;
#line 899 "zsytf2_rk.f"
			i__3 = j + k * a_dim1;
#line 899 "zsytf2_rk.f"
			z__2.r = z__3.r - a[i__3].r, z__2.i = z__3.i - a[i__3]
				.i;
#line 899 "zsytf2_rk.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 899 "zsytf2_rk.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;

/*                    Perform a rank-2 update of A(k+2:n,k+2:n) */

#line 903 "zsytf2_rk.f"
			i__2 = *n;
#line 903 "zsytf2_rk.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 904 "zsytf2_rk.f"
			    i__3 = i__ + j * a_dim1;
#line 904 "zsytf2_rk.f"
			    i__4 = i__ + j * a_dim1;
#line 904 "zsytf2_rk.f"
			    z_div(&z__4, &a[i__ + k * a_dim1], &d21);
#line 904 "zsytf2_rk.f"
			    z__3.r = z__4.r * wk.r - z__4.i * wk.i, z__3.i = 
				    z__4.r * wk.i + z__4.i * wk.r;
#line 904 "zsytf2_rk.f"
			    z__2.r = a[i__4].r - z__3.r, z__2.i = a[i__4].i - 
				    z__3.i;
#line 904 "zsytf2_rk.f"
			    z_div(&z__6, &a[i__ + (k + 1) * a_dim1], &d21);
#line 904 "zsytf2_rk.f"
			    z__5.r = z__6.r * wkp1.r - z__6.i * wkp1.i, 
				    z__5.i = z__6.r * wkp1.i + z__6.i * 
				    wkp1.r;
#line 904 "zsytf2_rk.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 904 "zsytf2_rk.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 906 "zsytf2_rk.f"
/* L50: */
#line 906 "zsytf2_rk.f"
			}

/*                    Store L(k) and L(k+1) in cols k and k+1 for row J */

#line 910 "zsytf2_rk.f"
			i__2 = j + k * a_dim1;
#line 910 "zsytf2_rk.f"
			z_div(&z__1, &wk, &d21);
#line 910 "zsytf2_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 911 "zsytf2_rk.f"
			i__2 = j + (k + 1) * a_dim1;
#line 911 "zsytf2_rk.f"
			z_div(&z__1, &wkp1, &d21);
#line 911 "zsytf2_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;

#line 913 "zsytf2_rk.f"
/* L60: */
#line 913 "zsytf2_rk.f"
		    }

#line 915 "zsytf2_rk.f"
		}

/*              Copy subdiagonal elements of D(K) to E(K) and */
/*              ZERO out subdiagonal entry of A */

#line 920 "zsytf2_rk.f"
		i__1 = k;
#line 920 "zsytf2_rk.f"
		i__2 = k + 1 + k * a_dim1;
#line 920 "zsytf2_rk.f"
		e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 921 "zsytf2_rk.f"
		i__1 = k + 1;
#line 921 "zsytf2_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 922 "zsytf2_rk.f"
		i__1 = k + 1 + k * a_dim1;
#line 922 "zsytf2_rk.f"
		a[i__1].r = 0., a[i__1].i = 0.;

#line 924 "zsytf2_rk.f"
	    }

/*           End column K is nonsingular */

#line 928 "zsytf2_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 932 "zsytf2_rk.f"
	if (kstep == 1) {
#line 933 "zsytf2_rk.f"
	    ipiv[k] = kp;
#line 934 "zsytf2_rk.f"
	} else {
#line 935 "zsytf2_rk.f"
	    ipiv[k] = -p;
#line 936 "zsytf2_rk.f"
	    ipiv[k + 1] = -kp;
#line 937 "zsytf2_rk.f"
	}

/*        Increase K and return to the start of the main loop */

#line 941 "zsytf2_rk.f"
	k += kstep;
#line 942 "zsytf2_rk.f"
	goto L40;

#line 944 "zsytf2_rk.f"
L64:

#line 946 "zsytf2_rk.f"
	;
#line 946 "zsytf2_rk.f"
    }

#line 948 "zsytf2_rk.f"
    return 0;

/*     End of ZSYTF2_RK */

} /* zsytf2_rk__ */

