#line 1 "zhetf2_rk.f"
/* zhetf2_rk.f -- translated by f2c (version 20100827).
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

#line 1 "zhetf2_rk.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZHETF2_RK computes the factorization of a complex Hermitian indefinite matrix using the bounded
 Bunch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETF2_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetf2_
rk.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetf2_
rk.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetf2_
rk.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETF2_RK( UPLO, N, A, LDA, E, IPIV, INFO ) */

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
/* > ZHETF2_RK computes the factorization of a complex Hermitian matrix A */
/* > using the bounded Bunch-Kaufman (rook) diagonal pivoting method: */
/* > */
/* >    A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is Hermitian and block */
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
/* >          Hermitian matrix A is stored: */
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
/* >          On entry, the Hermitian matrix A. */
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
/* >            a) ONLY diagonal elements of the Hermitian block diagonal */
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
/* >          elements of the Hermitian block diagonal matrix D */
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
/* >          the Hermitian block diagonal matrix D with 1-by-1 or 2-by-2 */
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

/* > \ingroup complex16HEcomputational */

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
/* Subroutine */ int zhetf2_rk__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *e, integer *ipiv, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k, p;
    static doublecomplex t;
    static doublereal r1, d11;
    static doublecomplex d12;
    static doublereal d22;
    static doublecomplex d21;
    static integer ii, kk, kp;
    static doublecomplex wk;
    static doublereal tt;
    static doublecomplex wkm1, wkp1;
    static logical done;
    static integer imax, jmax;
    extern /* Subroutine */ int zher_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dtemp, sfmin;
    static integer itemp, kstep;
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    static doublereal absakk;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
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

/*  ====================================================================== */

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

#line 299 "zhetf2_rk.f"
    /* Parameter adjustments */
#line 299 "zhetf2_rk.f"
    a_dim1 = *lda;
#line 299 "zhetf2_rk.f"
    a_offset = 1 + a_dim1;
#line 299 "zhetf2_rk.f"
    a -= a_offset;
#line 299 "zhetf2_rk.f"
    --e;
#line 299 "zhetf2_rk.f"
    --ipiv;
#line 299 "zhetf2_rk.f"

#line 299 "zhetf2_rk.f"
    /* Function Body */
#line 299 "zhetf2_rk.f"
    *info = 0;
#line 300 "zhetf2_rk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 301 "zhetf2_rk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 302 "zhetf2_rk.f"
	*info = -1;
#line 303 "zhetf2_rk.f"
    } else if (*n < 0) {
#line 304 "zhetf2_rk.f"
	*info = -2;
#line 305 "zhetf2_rk.f"
    } else if (*lda < max(1,*n)) {
#line 306 "zhetf2_rk.f"
	*info = -4;
#line 307 "zhetf2_rk.f"
    }
#line 308 "zhetf2_rk.f"
    if (*info != 0) {
#line 309 "zhetf2_rk.f"
	i__1 = -(*info);
#line 309 "zhetf2_rk.f"
	xerbla_("ZHETF2_RK", &i__1, (ftnlen)9);
#line 310 "zhetf2_rk.f"
	return 0;
#line 311 "zhetf2_rk.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 315 "zhetf2_rk.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 319 "zhetf2_rk.f"
    sfmin = dlamch_("S", (ftnlen)1);

#line 321 "zhetf2_rk.f"
    if (upper) {

/*        Factorize A as U*D*U**H using the upper triangle of A */

/*        Initilize the first entry of array E, where superdiagonal */
/*        elements of D are stored */

#line 328 "zhetf2_rk.f"
	e[1].r = 0., e[1].i = 0.;

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 333 "zhetf2_rk.f"
	k = *n;
#line 334 "zhetf2_rk.f"
L10:

/*        If K < 1, exit from loop */

#line 338 "zhetf2_rk.f"
	if (k < 1) {
#line 338 "zhetf2_rk.f"
	    goto L34;
#line 338 "zhetf2_rk.f"
	}
#line 340 "zhetf2_rk.f"
	kstep = 1;
#line 341 "zhetf2_rk.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 346 "zhetf2_rk.f"
	i__1 = k + k * a_dim1;
#line 346 "zhetf2_rk.f"
	absakk = (d__1 = a[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 352 "zhetf2_rk.f"
	if (k > 1) {
#line 353 "zhetf2_rk.f"
	    i__1 = k - 1;
#line 353 "zhetf2_rk.f"
	    imax = izamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 354 "zhetf2_rk.f"
	    i__1 = imax + k * a_dim1;
#line 354 "zhetf2_rk.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 355 "zhetf2_rk.f"
	} else {
#line 356 "zhetf2_rk.f"
	    colmax = 0.;
#line 357 "zhetf2_rk.f"
	}

#line 359 "zhetf2_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 363 "zhetf2_rk.f"
	    if (*info == 0) {
#line 363 "zhetf2_rk.f"
		*info = k;
#line 363 "zhetf2_rk.f"
	    }
#line 365 "zhetf2_rk.f"
	    kp = k;
#line 366 "zhetf2_rk.f"
	    i__1 = k + k * a_dim1;
#line 366 "zhetf2_rk.f"
	    i__2 = k + k * a_dim1;
#line 366 "zhetf2_rk.f"
	    d__1 = a[i__2].r;
#line 366 "zhetf2_rk.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;

/*           Set E( K ) to zero */

#line 370 "zhetf2_rk.f"
	    if (k > 1) {
#line 370 "zhetf2_rk.f"
		i__1 = k;
#line 370 "zhetf2_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 370 "zhetf2_rk.f"
	    }

#line 373 "zhetf2_rk.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 383 "zhetf2_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 387 "zhetf2_rk.f"
		kp = k;

#line 389 "zhetf2_rk.f"
	    } else {

#line 391 "zhetf2_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 395 "zhetf2_rk.f"
L12:

/*                 BEGIN pivot search loop body */


/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 404 "zhetf2_rk.f"
		if (imax != k) {
#line 405 "zhetf2_rk.f"
		    i__1 = k - imax;
#line 405 "zhetf2_rk.f"
		    jmax = imax + izamax_(&i__1, &a[imax + (imax + 1) * 
			    a_dim1], lda);
#line 407 "zhetf2_rk.f"
		    i__1 = imax + jmax * a_dim1;
#line 407 "zhetf2_rk.f"
		    rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    a[imax + jmax * a_dim1]), abs(d__2));
#line 408 "zhetf2_rk.f"
		} else {
#line 409 "zhetf2_rk.f"
		    rowmax = 0.;
#line 410 "zhetf2_rk.f"
		}

#line 412 "zhetf2_rk.f"
		if (imax > 1) {
#line 413 "zhetf2_rk.f"
		    i__1 = imax - 1;
#line 413 "zhetf2_rk.f"
		    itemp = izamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
#line 414 "zhetf2_rk.f"
		    i__1 = itemp + imax * a_dim1;
#line 414 "zhetf2_rk.f"
		    dtemp = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    itemp + imax * a_dim1]), abs(d__2));
#line 415 "zhetf2_rk.f"
		    if (dtemp > rowmax) {
#line 416 "zhetf2_rk.f"
			rowmax = dtemp;
#line 417 "zhetf2_rk.f"
			jmax = itemp;
#line 418 "zhetf2_rk.f"
		    }
#line 419 "zhetf2_rk.f"
		}

/*                 Case(2) */
/*                 Equivalent to testing for */
/*                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 426 "zhetf2_rk.f"
		i__1 = imax + imax * a_dim1;
#line 426 "zhetf2_rk.f"
		if (! ((d__1 = a[i__1].r, abs(d__1)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 432 "zhetf2_rk.f"
		    kp = imax;
#line 433 "zhetf2_rk.f"
		    done = TRUE_;

/*                 Case(3) */
/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 439 "zhetf2_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 445 "zhetf2_rk.f"
		    kp = imax;
#line 446 "zhetf2_rk.f"
		    kstep = 2;
#line 447 "zhetf2_rk.f"
		    done = TRUE_;

/*                 Case(4) */
#line 450 "zhetf2_rk.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 454 "zhetf2_rk.f"
		    p = imax;
#line 455 "zhetf2_rk.f"
		    colmax = rowmax;
#line 456 "zhetf2_rk.f"
		    imax = jmax;
#line 457 "zhetf2_rk.f"
		}

/*                 END pivot search loop body */

#line 461 "zhetf2_rk.f"
		if (! done) {
#line 461 "zhetf2_rk.f"
		    goto L12;
#line 461 "zhetf2_rk.f"
		}

#line 463 "zhetf2_rk.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 471 "zhetf2_rk.f"
	    kk = k - kstep + 1;

/*           For only a 2x2 pivot, interchange rows and columns K and P */
/*           in the leading submatrix A(1:k,1:k) */

#line 476 "zhetf2_rk.f"
	    if (kstep == 2 && p != k) {
/*              (1) Swap columnar parts */
#line 478 "zhetf2_rk.f"
		if (p > 1) {
#line 478 "zhetf2_rk.f"
		    i__1 = p - 1;
#line 478 "zhetf2_rk.f"
		    zswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 
			    1], &c__1);
#line 478 "zhetf2_rk.f"
		}
/*              (2) Swap and conjugate middle parts */
#line 481 "zhetf2_rk.f"
		i__1 = k - 1;
#line 481 "zhetf2_rk.f"
		for (j = p + 1; j <= i__1; ++j) {
#line 482 "zhetf2_rk.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 482 "zhetf2_rk.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 483 "zhetf2_rk.f"
		    i__2 = j + k * a_dim1;
#line 483 "zhetf2_rk.f"
		    d_cnjg(&z__1, &a[p + j * a_dim1]);
#line 483 "zhetf2_rk.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 484 "zhetf2_rk.f"
		    i__2 = p + j * a_dim1;
#line 484 "zhetf2_rk.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 485 "zhetf2_rk.f"
/* L14: */
#line 485 "zhetf2_rk.f"
		}
/*              (3) Swap and conjugate corner elements at row-col interserction */
#line 487 "zhetf2_rk.f"
		i__1 = p + k * a_dim1;
#line 487 "zhetf2_rk.f"
		d_cnjg(&z__1, &a[p + k * a_dim1]);
#line 487 "zhetf2_rk.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*              (4) Swap diagonal elements at row-col intersection */
#line 489 "zhetf2_rk.f"
		i__1 = k + k * a_dim1;
#line 489 "zhetf2_rk.f"
		r1 = a[i__1].r;
#line 490 "zhetf2_rk.f"
		i__1 = k + k * a_dim1;
#line 490 "zhetf2_rk.f"
		i__2 = p + p * a_dim1;
#line 490 "zhetf2_rk.f"
		d__1 = a[i__2].r;
#line 490 "zhetf2_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 491 "zhetf2_rk.f"
		i__1 = p + p * a_dim1;
#line 491 "zhetf2_rk.f"
		a[i__1].r = r1, a[i__1].i = 0.;

/*              Convert upper triangle of A into U form by applying */
/*              the interchanges in columns k+1:N. */

#line 496 "zhetf2_rk.f"
		if (k < *n) {
#line 496 "zhetf2_rk.f"
		    i__1 = *n - k;
#line 496 "zhetf2_rk.f"
		    zswap_(&i__1, &a[k + (k + 1) * a_dim1], lda, &a[p + (k + 
			    1) * a_dim1], lda);
#line 496 "zhetf2_rk.f"
		}

#line 499 "zhetf2_rk.f"
	    }

/*           For both 1x1 and 2x2 pivots, interchange rows and */
/*           columns KK and KP in the leading submatrix A(1:k,1:k) */

#line 504 "zhetf2_rk.f"
	    if (kp != kk) {
/*              (1) Swap columnar parts */
#line 506 "zhetf2_rk.f"
		if (kp > 1) {
#line 506 "zhetf2_rk.f"
		    i__1 = kp - 1;
#line 506 "zhetf2_rk.f"
		    zswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 506 "zhetf2_rk.f"
		}
/*              (2) Swap and conjugate middle parts */
#line 509 "zhetf2_rk.f"
		i__1 = kk - 1;
#line 509 "zhetf2_rk.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 510 "zhetf2_rk.f"
		    d_cnjg(&z__1, &a[j + kk * a_dim1]);
#line 510 "zhetf2_rk.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 511 "zhetf2_rk.f"
		    i__2 = j + kk * a_dim1;
#line 511 "zhetf2_rk.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 511 "zhetf2_rk.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 512 "zhetf2_rk.f"
		    i__2 = kp + j * a_dim1;
#line 512 "zhetf2_rk.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 513 "zhetf2_rk.f"
/* L15: */
#line 513 "zhetf2_rk.f"
		}
/*              (3) Swap and conjugate corner elements at row-col interserction */
#line 515 "zhetf2_rk.f"
		i__1 = kp + kk * a_dim1;
#line 515 "zhetf2_rk.f"
		d_cnjg(&z__1, &a[kp + kk * a_dim1]);
#line 515 "zhetf2_rk.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*              (4) Swap diagonal elements at row-col intersection */
#line 517 "zhetf2_rk.f"
		i__1 = kk + kk * a_dim1;
#line 517 "zhetf2_rk.f"
		r1 = a[i__1].r;
#line 518 "zhetf2_rk.f"
		i__1 = kk + kk * a_dim1;
#line 518 "zhetf2_rk.f"
		i__2 = kp + kp * a_dim1;
#line 518 "zhetf2_rk.f"
		d__1 = a[i__2].r;
#line 518 "zhetf2_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 519 "zhetf2_rk.f"
		i__1 = kp + kp * a_dim1;
#line 519 "zhetf2_rk.f"
		a[i__1].r = r1, a[i__1].i = 0.;

#line 521 "zhetf2_rk.f"
		if (kstep == 2) {
/*                 (*) Make sure that diagonal element of pivot is real */
#line 523 "zhetf2_rk.f"
		    i__1 = k + k * a_dim1;
#line 523 "zhetf2_rk.f"
		    i__2 = k + k * a_dim1;
#line 523 "zhetf2_rk.f"
		    d__1 = a[i__2].r;
#line 523 "zhetf2_rk.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
/*                 (5) Swap row elements */
#line 525 "zhetf2_rk.f"
		    i__1 = k - 1 + k * a_dim1;
#line 525 "zhetf2_rk.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 526 "zhetf2_rk.f"
		    i__1 = k - 1 + k * a_dim1;
#line 526 "zhetf2_rk.f"
		    i__2 = kp + k * a_dim1;
#line 526 "zhetf2_rk.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 527 "zhetf2_rk.f"
		    i__1 = kp + k * a_dim1;
#line 527 "zhetf2_rk.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 528 "zhetf2_rk.f"
		}

/*              Convert upper triangle of A into U form by applying */
/*              the interchanges in columns k+1:N. */

#line 533 "zhetf2_rk.f"
		if (k < *n) {
#line 533 "zhetf2_rk.f"
		    i__1 = *n - k;
#line 533 "zhetf2_rk.f"
		    zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 533 "zhetf2_rk.f"
		}

#line 537 "zhetf2_rk.f"
	    } else {
/*              (*) Make sure that diagonal element of pivot is real */
#line 539 "zhetf2_rk.f"
		i__1 = k + k * a_dim1;
#line 539 "zhetf2_rk.f"
		i__2 = k + k * a_dim1;
#line 539 "zhetf2_rk.f"
		d__1 = a[i__2].r;
#line 539 "zhetf2_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 540 "zhetf2_rk.f"
		if (kstep == 2) {
#line 540 "zhetf2_rk.f"
		    i__1 = k - 1 + (k - 1) * a_dim1;
#line 540 "zhetf2_rk.f"
		    i__2 = k - 1 + (k - 1) * a_dim1;
#line 540 "zhetf2_rk.f"
		    d__1 = a[i__2].r;
#line 540 "zhetf2_rk.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 540 "zhetf2_rk.f"
		}
#line 542 "zhetf2_rk.f"
	    }

/*           Update the leading submatrix */

#line 546 "zhetf2_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

#line 554 "zhetf2_rk.f"
		if (k > 1) {

/*                 Perform a rank-1 update of A(1:k-1,1:k-1) and */
/*                 store U(k) in column k */

#line 559 "zhetf2_rk.f"
		    i__1 = k + k * a_dim1;
#line 559 "zhetf2_rk.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) >= sfmin) {

/*                    Perform a rank-1 update of A(1:k-1,1:k-1) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*1/D(k)*W(k)**T */

#line 565 "zhetf2_rk.f"
			i__1 = k + k * a_dim1;
#line 565 "zhetf2_rk.f"
			d11 = 1. / a[i__1].r;
#line 566 "zhetf2_rk.f"
			i__1 = k - 1;
#line 566 "zhetf2_rk.f"
			d__1 = -d11;
#line 566 "zhetf2_rk.f"
			zher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);

/*                    Store U(k) in column k */

#line 570 "zhetf2_rk.f"
			i__1 = k - 1;
#line 570 "zhetf2_rk.f"
			zdscal_(&i__1, &d11, &a[k * a_dim1 + 1], &c__1);
#line 571 "zhetf2_rk.f"
		    } else {

/*                    Store L(k) in column K */

#line 575 "zhetf2_rk.f"
			i__1 = k + k * a_dim1;
#line 575 "zhetf2_rk.f"
			d11 = a[i__1].r;
#line 576 "zhetf2_rk.f"
			i__1 = k - 1;
#line 576 "zhetf2_rk.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 577 "zhetf2_rk.f"
			    i__2 = ii + k * a_dim1;
#line 577 "zhetf2_rk.f"
			    i__3 = ii + k * a_dim1;
#line 577 "zhetf2_rk.f"
			    z__1.r = a[i__3].r / d11, z__1.i = a[i__3].i / 
				    d11;
#line 577 "zhetf2_rk.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 578 "zhetf2_rk.f"
/* L16: */
#line 578 "zhetf2_rk.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 585 "zhetf2_rk.f"
			i__1 = k - 1;
#line 585 "zhetf2_rk.f"
			d__1 = -d11;
#line 585 "zhetf2_rk.f"
			zher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);
#line 586 "zhetf2_rk.f"
		    }

/*                 Store the superdiagonal element of D in array E */

#line 590 "zhetf2_rk.f"
		    i__1 = k;
#line 590 "zhetf2_rk.f"
		    e[i__1].r = 0., e[i__1].i = 0.;

#line 592 "zhetf2_rk.f"
		}

#line 594 "zhetf2_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 610 "zhetf2_rk.f"
		if (k > 2) {
/*                 D = |A12| */
#line 612 "zhetf2_rk.f"
		    i__1 = k - 1 + k * a_dim1;
#line 612 "zhetf2_rk.f"
		    d__1 = a[i__1].r;
#line 612 "zhetf2_rk.f"
		    d__2 = d_imag(&a[k - 1 + k * a_dim1]);
#line 612 "zhetf2_rk.f"
		    d__ = dlapy2_(&d__1, &d__2);
#line 614 "zhetf2_rk.f"
		    i__1 = k + k * a_dim1;
#line 614 "zhetf2_rk.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 614 "zhetf2_rk.f"
		    d11 = z__1.r;
#line 615 "zhetf2_rk.f"
		    i__1 = k - 1 + (k - 1) * a_dim1;
#line 615 "zhetf2_rk.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 615 "zhetf2_rk.f"
		    d22 = z__1.r;
#line 616 "zhetf2_rk.f"
		    i__1 = k - 1 + k * a_dim1;
#line 616 "zhetf2_rk.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 616 "zhetf2_rk.f"
		    d12.r = z__1.r, d12.i = z__1.i;
#line 617 "zhetf2_rk.f"
		    tt = 1. / (d11 * d22 - 1.);

#line 619 "zhetf2_rk.f"
		    for (j = k - 2; j >= 1; --j) {

/*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */

#line 623 "zhetf2_rk.f"
			i__1 = j + (k - 1) * a_dim1;
#line 623 "zhetf2_rk.f"
			z__3.r = d11 * a[i__1].r, z__3.i = d11 * a[i__1].i;
#line 623 "zhetf2_rk.f"
			d_cnjg(&z__5, &d12);
#line 623 "zhetf2_rk.f"
			i__2 = j + k * a_dim1;
#line 623 "zhetf2_rk.f"
			z__4.r = z__5.r * a[i__2].r - z__5.i * a[i__2].i, 
				z__4.i = z__5.r * a[i__2].i + z__5.i * a[i__2]
				.r;
#line 623 "zhetf2_rk.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 623 "zhetf2_rk.f"
			z__1.r = tt * z__2.r, z__1.i = tt * z__2.i;
#line 623 "zhetf2_rk.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 625 "zhetf2_rk.f"
			i__1 = j + k * a_dim1;
#line 625 "zhetf2_rk.f"
			z__3.r = d22 * a[i__1].r, z__3.i = d22 * a[i__1].i;
#line 625 "zhetf2_rk.f"
			i__2 = j + (k - 1) * a_dim1;
#line 625 "zhetf2_rk.f"
			z__4.r = d12.r * a[i__2].r - d12.i * a[i__2].i, 
				z__4.i = d12.r * a[i__2].i + d12.i * a[i__2]
				.r;
#line 625 "zhetf2_rk.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 625 "zhetf2_rk.f"
			z__1.r = tt * z__2.r, z__1.i = tt * z__2.i;
#line 625 "zhetf2_rk.f"
			wk.r = z__1.r, wk.i = z__1.i;

/*                    Perform a rank-2 update of A(1:k-2,1:k-2) */

#line 629 "zhetf2_rk.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 630 "zhetf2_rk.f"
			    i__1 = i__ + j * a_dim1;
#line 630 "zhetf2_rk.f"
			    i__2 = i__ + j * a_dim1;
#line 630 "zhetf2_rk.f"
			    i__3 = i__ + k * a_dim1;
#line 630 "zhetf2_rk.f"
			    z__4.r = a[i__3].r / d__, z__4.i = a[i__3].i / 
				    d__;
#line 630 "zhetf2_rk.f"
			    d_cnjg(&z__5, &wk);
#line 630 "zhetf2_rk.f"
			    z__3.r = z__4.r * z__5.r - z__4.i * z__5.i, 
				    z__3.i = z__4.r * z__5.i + z__4.i * 
				    z__5.r;
#line 630 "zhetf2_rk.f"
			    z__2.r = a[i__2].r - z__3.r, z__2.i = a[i__2].i - 
				    z__3.i;
#line 630 "zhetf2_rk.f"
			    i__4 = i__ + (k - 1) * a_dim1;
#line 630 "zhetf2_rk.f"
			    z__7.r = a[i__4].r / d__, z__7.i = a[i__4].i / 
				    d__;
#line 630 "zhetf2_rk.f"
			    d_cnjg(&z__8, &wkm1);
#line 630 "zhetf2_rk.f"
			    z__6.r = z__7.r * z__8.r - z__7.i * z__8.i, 
				    z__6.i = z__7.r * z__8.i + z__7.i * 
				    z__8.r;
#line 630 "zhetf2_rk.f"
			    z__1.r = z__2.r - z__6.r, z__1.i = z__2.i - 
				    z__6.i;
#line 630 "zhetf2_rk.f"
			    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 633 "zhetf2_rk.f"
/* L20: */
#line 633 "zhetf2_rk.f"
			}

/*                    Store U(k) and U(k-1) in cols k and k-1 for row J */

#line 637 "zhetf2_rk.f"
			i__1 = j + k * a_dim1;
#line 637 "zhetf2_rk.f"
			z__1.r = wk.r / d__, z__1.i = wk.i / d__;
#line 637 "zhetf2_rk.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 638 "zhetf2_rk.f"
			i__1 = j + (k - 1) * a_dim1;
#line 638 "zhetf2_rk.f"
			z__1.r = wkm1.r / d__, z__1.i = wkm1.i / d__;
#line 638 "zhetf2_rk.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*                    (*) Make sure that diagonal element of pivot is real */
#line 640 "zhetf2_rk.f"
			i__1 = j + j * a_dim1;
#line 640 "zhetf2_rk.f"
			i__2 = j + j * a_dim1;
#line 640 "zhetf2_rk.f"
			d__1 = a[i__2].r;
#line 640 "zhetf2_rk.f"
			z__1.r = d__1, z__1.i = 0.;
#line 640 "zhetf2_rk.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 642 "zhetf2_rk.f"
/* L30: */
#line 642 "zhetf2_rk.f"
		    }

#line 644 "zhetf2_rk.f"
		}

/*              Copy superdiagonal elements of D(K) to E(K) and */
/*              ZERO out superdiagonal entry of A */

#line 649 "zhetf2_rk.f"
		i__1 = k;
#line 649 "zhetf2_rk.f"
		i__2 = k - 1 + k * a_dim1;
#line 649 "zhetf2_rk.f"
		e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 650 "zhetf2_rk.f"
		i__1 = k - 1;
#line 650 "zhetf2_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 651 "zhetf2_rk.f"
		i__1 = k - 1 + k * a_dim1;
#line 651 "zhetf2_rk.f"
		a[i__1].r = 0., a[i__1].i = 0.;

#line 653 "zhetf2_rk.f"
	    }

/*           End column K is nonsingular */

#line 657 "zhetf2_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 661 "zhetf2_rk.f"
	if (kstep == 1) {
#line 662 "zhetf2_rk.f"
	    ipiv[k] = kp;
#line 663 "zhetf2_rk.f"
	} else {
#line 664 "zhetf2_rk.f"
	    ipiv[k] = -p;
#line 665 "zhetf2_rk.f"
	    ipiv[k - 1] = -kp;
#line 666 "zhetf2_rk.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 670 "zhetf2_rk.f"
	k -= kstep;
#line 671 "zhetf2_rk.f"
	goto L10;

#line 673 "zhetf2_rk.f"
L34:

#line 675 "zhetf2_rk.f"
	;
#line 675 "zhetf2_rk.f"
    } else {

/*        Factorize A as L*D*L**H using the lower triangle of A */

/*        Initilize the unused last entry of the subdiagonal array E. */

#line 681 "zhetf2_rk.f"
	i__1 = *n;
#line 681 "zhetf2_rk.f"
	e[i__1].r = 0., e[i__1].i = 0.;

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 686 "zhetf2_rk.f"
	k = 1;
#line 687 "zhetf2_rk.f"
L40:

/*        If K > N, exit from loop */

#line 691 "zhetf2_rk.f"
	if (k > *n) {
#line 691 "zhetf2_rk.f"
	    goto L64;
#line 691 "zhetf2_rk.f"
	}
#line 693 "zhetf2_rk.f"
	kstep = 1;
#line 694 "zhetf2_rk.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 699 "zhetf2_rk.f"
	i__1 = k + k * a_dim1;
#line 699 "zhetf2_rk.f"
	absakk = (d__1 = a[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 705 "zhetf2_rk.f"
	if (k < *n) {
#line 706 "zhetf2_rk.f"
	    i__1 = *n - k;
#line 706 "zhetf2_rk.f"
	    imax = k + izamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 707 "zhetf2_rk.f"
	    i__1 = imax + k * a_dim1;
#line 707 "zhetf2_rk.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 708 "zhetf2_rk.f"
	} else {
#line 709 "zhetf2_rk.f"
	    colmax = 0.;
#line 710 "zhetf2_rk.f"
	}

#line 712 "zhetf2_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 716 "zhetf2_rk.f"
	    if (*info == 0) {
#line 716 "zhetf2_rk.f"
		*info = k;
#line 716 "zhetf2_rk.f"
	    }
#line 718 "zhetf2_rk.f"
	    kp = k;
#line 719 "zhetf2_rk.f"
	    i__1 = k + k * a_dim1;
#line 719 "zhetf2_rk.f"
	    i__2 = k + k * a_dim1;
#line 719 "zhetf2_rk.f"
	    d__1 = a[i__2].r;
#line 719 "zhetf2_rk.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;

/*           Set E( K ) to zero */

#line 723 "zhetf2_rk.f"
	    if (k < *n) {
#line 723 "zhetf2_rk.f"
		i__1 = k;
#line 723 "zhetf2_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 723 "zhetf2_rk.f"
	    }

#line 726 "zhetf2_rk.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 736 "zhetf2_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 740 "zhetf2_rk.f"
		kp = k;

#line 742 "zhetf2_rk.f"
	    } else {

#line 744 "zhetf2_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 748 "zhetf2_rk.f"
L42:

/*                 BEGIN pivot search loop body */


/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 757 "zhetf2_rk.f"
		if (imax != k) {
#line 758 "zhetf2_rk.f"
		    i__1 = imax - k;
#line 758 "zhetf2_rk.f"
		    jmax = k - 1 + izamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 759 "zhetf2_rk.f"
		    i__1 = imax + jmax * a_dim1;
#line 759 "zhetf2_rk.f"
		    rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    a[imax + jmax * a_dim1]), abs(d__2));
#line 760 "zhetf2_rk.f"
		} else {
#line 761 "zhetf2_rk.f"
		    rowmax = 0.;
#line 762 "zhetf2_rk.f"
		}

#line 764 "zhetf2_rk.f"
		if (imax < *n) {
#line 765 "zhetf2_rk.f"
		    i__1 = *n - imax;
#line 765 "zhetf2_rk.f"
		    itemp = imax + izamax_(&i__1, &a[imax + 1 + imax * a_dim1]
			    , &c__1);
#line 767 "zhetf2_rk.f"
		    i__1 = itemp + imax * a_dim1;
#line 767 "zhetf2_rk.f"
		    dtemp = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    itemp + imax * a_dim1]), abs(d__2));
#line 768 "zhetf2_rk.f"
		    if (dtemp > rowmax) {
#line 769 "zhetf2_rk.f"
			rowmax = dtemp;
#line 770 "zhetf2_rk.f"
			jmax = itemp;
#line 771 "zhetf2_rk.f"
		    }
#line 772 "zhetf2_rk.f"
		}

/*                 Case(2) */
/*                 Equivalent to testing for */
/*                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 779 "zhetf2_rk.f"
		i__1 = imax + imax * a_dim1;
#line 779 "zhetf2_rk.f"
		if (! ((d__1 = a[i__1].r, abs(d__1)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 785 "zhetf2_rk.f"
		    kp = imax;
#line 786 "zhetf2_rk.f"
		    done = TRUE_;

/*                 Case(3) */
/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 792 "zhetf2_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 798 "zhetf2_rk.f"
		    kp = imax;
#line 799 "zhetf2_rk.f"
		    kstep = 2;
#line 800 "zhetf2_rk.f"
		    done = TRUE_;

/*                 Case(4) */
#line 803 "zhetf2_rk.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 807 "zhetf2_rk.f"
		    p = imax;
#line 808 "zhetf2_rk.f"
		    colmax = rowmax;
#line 809 "zhetf2_rk.f"
		    imax = jmax;
#line 810 "zhetf2_rk.f"
		}


/*                 END pivot search loop body */

#line 815 "zhetf2_rk.f"
		if (! done) {
#line 815 "zhetf2_rk.f"
		    goto L42;
#line 815 "zhetf2_rk.f"
		}

#line 817 "zhetf2_rk.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 825 "zhetf2_rk.f"
	    kk = k + kstep - 1;

/*           For only a 2x2 pivot, interchange rows and columns K and P */
/*           in the trailing submatrix A(k:n,k:n) */

#line 830 "zhetf2_rk.f"
	    if (kstep == 2 && p != k) {
/*              (1) Swap columnar parts */
#line 832 "zhetf2_rk.f"
		if (p < *n) {
#line 832 "zhetf2_rk.f"
		    i__1 = *n - p;
#line 832 "zhetf2_rk.f"
		    zswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p 
			    * a_dim1], &c__1);
#line 832 "zhetf2_rk.f"
		}
/*              (2) Swap and conjugate middle parts */
#line 835 "zhetf2_rk.f"
		i__1 = p - 1;
#line 835 "zhetf2_rk.f"
		for (j = k + 1; j <= i__1; ++j) {
#line 836 "zhetf2_rk.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 836 "zhetf2_rk.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 837 "zhetf2_rk.f"
		    i__2 = j + k * a_dim1;
#line 837 "zhetf2_rk.f"
		    d_cnjg(&z__1, &a[p + j * a_dim1]);
#line 837 "zhetf2_rk.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 838 "zhetf2_rk.f"
		    i__2 = p + j * a_dim1;
#line 838 "zhetf2_rk.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 839 "zhetf2_rk.f"
/* L44: */
#line 839 "zhetf2_rk.f"
		}
/*              (3) Swap and conjugate corner elements at row-col interserction */
#line 841 "zhetf2_rk.f"
		i__1 = p + k * a_dim1;
#line 841 "zhetf2_rk.f"
		d_cnjg(&z__1, &a[p + k * a_dim1]);
#line 841 "zhetf2_rk.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*              (4) Swap diagonal elements at row-col intersection */
#line 843 "zhetf2_rk.f"
		i__1 = k + k * a_dim1;
#line 843 "zhetf2_rk.f"
		r1 = a[i__1].r;
#line 844 "zhetf2_rk.f"
		i__1 = k + k * a_dim1;
#line 844 "zhetf2_rk.f"
		i__2 = p + p * a_dim1;
#line 844 "zhetf2_rk.f"
		d__1 = a[i__2].r;
#line 844 "zhetf2_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 845 "zhetf2_rk.f"
		i__1 = p + p * a_dim1;
#line 845 "zhetf2_rk.f"
		a[i__1].r = r1, a[i__1].i = 0.;

/*              Convert lower triangle of A into L form by applying */
/*              the interchanges in columns 1:k-1. */

#line 850 "zhetf2_rk.f"
		if (k > 1) {
#line 850 "zhetf2_rk.f"
		    i__1 = k - 1;
#line 850 "zhetf2_rk.f"
		    zswap_(&i__1, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 850 "zhetf2_rk.f"
		}

#line 853 "zhetf2_rk.f"
	    }

/*           For both 1x1 and 2x2 pivots, interchange rows and */
/*           columns KK and KP in the trailing submatrix A(k:n,k:n) */

#line 858 "zhetf2_rk.f"
	    if (kp != kk) {
/*              (1) Swap columnar parts */
#line 860 "zhetf2_rk.f"
		if (kp < *n) {
#line 860 "zhetf2_rk.f"
		    i__1 = *n - kp;
#line 860 "zhetf2_rk.f"
		    zswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 860 "zhetf2_rk.f"
		}
/*              (2) Swap and conjugate middle parts */
#line 863 "zhetf2_rk.f"
		i__1 = kp - 1;
#line 863 "zhetf2_rk.f"
		for (j = kk + 1; j <= i__1; ++j) {
#line 864 "zhetf2_rk.f"
		    d_cnjg(&z__1, &a[j + kk * a_dim1]);
#line 864 "zhetf2_rk.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 865 "zhetf2_rk.f"
		    i__2 = j + kk * a_dim1;
#line 865 "zhetf2_rk.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 865 "zhetf2_rk.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 866 "zhetf2_rk.f"
		    i__2 = kp + j * a_dim1;
#line 866 "zhetf2_rk.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 867 "zhetf2_rk.f"
/* L45: */
#line 867 "zhetf2_rk.f"
		}
/*              (3) Swap and conjugate corner elements at row-col interserction */
#line 869 "zhetf2_rk.f"
		i__1 = kp + kk * a_dim1;
#line 869 "zhetf2_rk.f"
		d_cnjg(&z__1, &a[kp + kk * a_dim1]);
#line 869 "zhetf2_rk.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*              (4) Swap diagonal elements at row-col intersection */
#line 871 "zhetf2_rk.f"
		i__1 = kk + kk * a_dim1;
#line 871 "zhetf2_rk.f"
		r1 = a[i__1].r;
#line 872 "zhetf2_rk.f"
		i__1 = kk + kk * a_dim1;
#line 872 "zhetf2_rk.f"
		i__2 = kp + kp * a_dim1;
#line 872 "zhetf2_rk.f"
		d__1 = a[i__2].r;
#line 872 "zhetf2_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 873 "zhetf2_rk.f"
		i__1 = kp + kp * a_dim1;
#line 873 "zhetf2_rk.f"
		a[i__1].r = r1, a[i__1].i = 0.;

#line 875 "zhetf2_rk.f"
		if (kstep == 2) {
/*                 (*) Make sure that diagonal element of pivot is real */
#line 877 "zhetf2_rk.f"
		    i__1 = k + k * a_dim1;
#line 877 "zhetf2_rk.f"
		    i__2 = k + k * a_dim1;
#line 877 "zhetf2_rk.f"
		    d__1 = a[i__2].r;
#line 877 "zhetf2_rk.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
/*                 (5) Swap row elements */
#line 879 "zhetf2_rk.f"
		    i__1 = k + 1 + k * a_dim1;
#line 879 "zhetf2_rk.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 880 "zhetf2_rk.f"
		    i__1 = k + 1 + k * a_dim1;
#line 880 "zhetf2_rk.f"
		    i__2 = kp + k * a_dim1;
#line 880 "zhetf2_rk.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 881 "zhetf2_rk.f"
		    i__1 = kp + k * a_dim1;
#line 881 "zhetf2_rk.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 882 "zhetf2_rk.f"
		}

/*              Convert lower triangle of A into L form by applying */
/*              the interchanges in columns 1:k-1. */

#line 887 "zhetf2_rk.f"
		if (k > 1) {
#line 887 "zhetf2_rk.f"
		    i__1 = k - 1;
#line 887 "zhetf2_rk.f"
		    zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 887 "zhetf2_rk.f"
		}

#line 890 "zhetf2_rk.f"
	    } else {
/*              (*) Make sure that diagonal element of pivot is real */
#line 892 "zhetf2_rk.f"
		i__1 = k + k * a_dim1;
#line 892 "zhetf2_rk.f"
		i__2 = k + k * a_dim1;
#line 892 "zhetf2_rk.f"
		d__1 = a[i__2].r;
#line 892 "zhetf2_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 893 "zhetf2_rk.f"
		if (kstep == 2) {
#line 893 "zhetf2_rk.f"
		    i__1 = k + 1 + (k + 1) * a_dim1;
#line 893 "zhetf2_rk.f"
		    i__2 = k + 1 + (k + 1) * a_dim1;
#line 893 "zhetf2_rk.f"
		    d__1 = a[i__2].r;
#line 893 "zhetf2_rk.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 893 "zhetf2_rk.f"
		}
#line 895 "zhetf2_rk.f"
	    }

/*           Update the trailing submatrix */

#line 899 "zhetf2_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of A now holds */

/*              W(k) = L(k)*D(k), */

/*              where L(k) is the k-th column of L */

#line 907 "zhetf2_rk.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) and */
/*                 store L(k) in column k */

/*                 Handle division by a small number */

#line 914 "zhetf2_rk.f"
		    i__1 = k + k * a_dim1;
#line 914 "zhetf2_rk.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) >= sfmin) {

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */

#line 920 "zhetf2_rk.f"
			i__1 = k + k * a_dim1;
#line 920 "zhetf2_rk.f"
			d11 = 1. / a[i__1].r;
#line 921 "zhetf2_rk.f"
			i__1 = *n - k;
#line 921 "zhetf2_rk.f"
			d__1 = -d11;
#line 921 "zhetf2_rk.f"
			zher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);

/*                    Store L(k) in column k */

#line 926 "zhetf2_rk.f"
			i__1 = *n - k;
#line 926 "zhetf2_rk.f"
			zdscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
#line 927 "zhetf2_rk.f"
		    } else {

/*                    Store L(k) in column k */

#line 931 "zhetf2_rk.f"
			i__1 = k + k * a_dim1;
#line 931 "zhetf2_rk.f"
			d11 = a[i__1].r;
#line 932 "zhetf2_rk.f"
			i__1 = *n;
#line 932 "zhetf2_rk.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 933 "zhetf2_rk.f"
			    i__2 = ii + k * a_dim1;
#line 933 "zhetf2_rk.f"
			    i__3 = ii + k * a_dim1;
#line 933 "zhetf2_rk.f"
			    z__1.r = a[i__3].r / d11, z__1.i = a[i__3].i / 
				    d11;
#line 933 "zhetf2_rk.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 934 "zhetf2_rk.f"
/* L46: */
#line 934 "zhetf2_rk.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 941 "zhetf2_rk.f"
			i__1 = *n - k;
#line 941 "zhetf2_rk.f"
			d__1 = -d11;
#line 941 "zhetf2_rk.f"
			zher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);
#line 943 "zhetf2_rk.f"
		    }

/*                 Store the subdiagonal element of D in array E */

#line 947 "zhetf2_rk.f"
		    i__1 = k;
#line 947 "zhetf2_rk.f"
		    e[i__1].r = 0., e[i__1].i = 0.;

#line 949 "zhetf2_rk.f"
		}

#line 951 "zhetf2_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */


/*              Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T */
/*                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 968 "zhetf2_rk.f"
		if (k < *n - 1) {
/*                 D = |A21| */
#line 970 "zhetf2_rk.f"
		    i__1 = k + 1 + k * a_dim1;
#line 970 "zhetf2_rk.f"
		    d__1 = a[i__1].r;
#line 970 "zhetf2_rk.f"
		    d__2 = d_imag(&a[k + 1 + k * a_dim1]);
#line 970 "zhetf2_rk.f"
		    d__ = dlapy2_(&d__1, &d__2);
#line 972 "zhetf2_rk.f"
		    i__1 = k + 1 + (k + 1) * a_dim1;
#line 972 "zhetf2_rk.f"
		    d11 = a[i__1].r / d__;
#line 973 "zhetf2_rk.f"
		    i__1 = k + k * a_dim1;
#line 973 "zhetf2_rk.f"
		    d22 = a[i__1].r / d__;
#line 974 "zhetf2_rk.f"
		    i__1 = k + 1 + k * a_dim1;
#line 974 "zhetf2_rk.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 974 "zhetf2_rk.f"
		    d21.r = z__1.r, d21.i = z__1.i;
#line 975 "zhetf2_rk.f"
		    tt = 1. / (d11 * d22 - 1.);

#line 977 "zhetf2_rk.f"
		    i__1 = *n;
#line 977 "zhetf2_rk.f"
		    for (j = k + 2; j <= i__1; ++j) {

/*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */

#line 981 "zhetf2_rk.f"
			i__2 = j + k * a_dim1;
#line 981 "zhetf2_rk.f"
			z__3.r = d11 * a[i__2].r, z__3.i = d11 * a[i__2].i;
#line 981 "zhetf2_rk.f"
			i__3 = j + (k + 1) * a_dim1;
#line 981 "zhetf2_rk.f"
			z__4.r = d21.r * a[i__3].r - d21.i * a[i__3].i, 
				z__4.i = d21.r * a[i__3].i + d21.i * a[i__3]
				.r;
#line 981 "zhetf2_rk.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 981 "zhetf2_rk.f"
			z__1.r = tt * z__2.r, z__1.i = tt * z__2.i;
#line 981 "zhetf2_rk.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 982 "zhetf2_rk.f"
			i__2 = j + (k + 1) * a_dim1;
#line 982 "zhetf2_rk.f"
			z__3.r = d22 * a[i__2].r, z__3.i = d22 * a[i__2].i;
#line 982 "zhetf2_rk.f"
			d_cnjg(&z__5, &d21);
#line 982 "zhetf2_rk.f"
			i__3 = j + k * a_dim1;
#line 982 "zhetf2_rk.f"
			z__4.r = z__5.r * a[i__3].r - z__5.i * a[i__3].i, 
				z__4.i = z__5.r * a[i__3].i + z__5.i * a[i__3]
				.r;
#line 982 "zhetf2_rk.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 982 "zhetf2_rk.f"
			z__1.r = tt * z__2.r, z__1.i = tt * z__2.i;
#line 982 "zhetf2_rk.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;

/*                    Perform a rank-2 update of A(k+2:n,k+2:n) */

#line 987 "zhetf2_rk.f"
			i__2 = *n;
#line 987 "zhetf2_rk.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 988 "zhetf2_rk.f"
			    i__3 = i__ + j * a_dim1;
#line 988 "zhetf2_rk.f"
			    i__4 = i__ + j * a_dim1;
#line 988 "zhetf2_rk.f"
			    i__5 = i__ + k * a_dim1;
#line 988 "zhetf2_rk.f"
			    z__4.r = a[i__5].r / d__, z__4.i = a[i__5].i / 
				    d__;
#line 988 "zhetf2_rk.f"
			    d_cnjg(&z__5, &wk);
#line 988 "zhetf2_rk.f"
			    z__3.r = z__4.r * z__5.r - z__4.i * z__5.i, 
				    z__3.i = z__4.r * z__5.i + z__4.i * 
				    z__5.r;
#line 988 "zhetf2_rk.f"
			    z__2.r = a[i__4].r - z__3.r, z__2.i = a[i__4].i - 
				    z__3.i;
#line 988 "zhetf2_rk.f"
			    i__6 = i__ + (k + 1) * a_dim1;
#line 988 "zhetf2_rk.f"
			    z__7.r = a[i__6].r / d__, z__7.i = a[i__6].i / 
				    d__;
#line 988 "zhetf2_rk.f"
			    d_cnjg(&z__8, &wkp1);
#line 988 "zhetf2_rk.f"
			    z__6.r = z__7.r * z__8.r - z__7.i * z__8.i, 
				    z__6.i = z__7.r * z__8.i + z__7.i * 
				    z__8.r;
#line 988 "zhetf2_rk.f"
			    z__1.r = z__2.r - z__6.r, z__1.i = z__2.i - 
				    z__6.i;
#line 988 "zhetf2_rk.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 991 "zhetf2_rk.f"
/* L50: */
#line 991 "zhetf2_rk.f"
			}

/*                    Store L(k) and L(k+1) in cols k and k+1 for row J */

#line 995 "zhetf2_rk.f"
			i__2 = j + k * a_dim1;
#line 995 "zhetf2_rk.f"
			z__1.r = wk.r / d__, z__1.i = wk.i / d__;
#line 995 "zhetf2_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 996 "zhetf2_rk.f"
			i__2 = j + (k + 1) * a_dim1;
#line 996 "zhetf2_rk.f"
			z__1.r = wkp1.r / d__, z__1.i = wkp1.i / d__;
#line 996 "zhetf2_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/*                    (*) Make sure that diagonal element of pivot is real */
#line 998 "zhetf2_rk.f"
			i__2 = j + j * a_dim1;
#line 998 "zhetf2_rk.f"
			i__3 = j + j * a_dim1;
#line 998 "zhetf2_rk.f"
			d__1 = a[i__3].r;
#line 998 "zhetf2_rk.f"
			z__1.r = d__1, z__1.i = 0.;
#line 998 "zhetf2_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;

#line 1000 "zhetf2_rk.f"
/* L60: */
#line 1000 "zhetf2_rk.f"
		    }

#line 1002 "zhetf2_rk.f"
		}

/*              Copy subdiagonal elements of D(K) to E(K) and */
/*              ZERO out subdiagonal entry of A */

#line 1007 "zhetf2_rk.f"
		i__1 = k;
#line 1007 "zhetf2_rk.f"
		i__2 = k + 1 + k * a_dim1;
#line 1007 "zhetf2_rk.f"
		e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 1008 "zhetf2_rk.f"
		i__1 = k + 1;
#line 1008 "zhetf2_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 1009 "zhetf2_rk.f"
		i__1 = k + 1 + k * a_dim1;
#line 1009 "zhetf2_rk.f"
		a[i__1].r = 0., a[i__1].i = 0.;

#line 1011 "zhetf2_rk.f"
	    }

/*           End column K is nonsingular */

#line 1015 "zhetf2_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 1019 "zhetf2_rk.f"
	if (kstep == 1) {
#line 1020 "zhetf2_rk.f"
	    ipiv[k] = kp;
#line 1021 "zhetf2_rk.f"
	} else {
#line 1022 "zhetf2_rk.f"
	    ipiv[k] = -p;
#line 1023 "zhetf2_rk.f"
	    ipiv[k + 1] = -kp;
#line 1024 "zhetf2_rk.f"
	}

/*        Increase K and return to the start of the main loop */

#line 1028 "zhetf2_rk.f"
	k += kstep;
#line 1029 "zhetf2_rk.f"
	goto L40;

#line 1031 "zhetf2_rk.f"
L64:

#line 1033 "zhetf2_rk.f"
	;
#line 1033 "zhetf2_rk.f"
    }

#line 1035 "zhetf2_rk.f"
    return 0;

/*     End of ZHETF2_RK */

} /* zhetf2_rk__ */

