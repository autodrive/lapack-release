#line 1 "slasyf_rk.f"
/* slasyf_rk.f -- translated by f2c (version 20100827).
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

#line 1 "slasyf_rk.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = -1.;
static doublereal c_b10 = 1.;

/* > \brief \b SLASYF_RK computes a partial factorization of a real symmetric indefinite matrix using bounded 
Bunch-Kaufman (rook) diagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASYF_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf_
rk.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf_
rk.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf_
rk.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASYF_RK( UPLO, N, NB, KB, A, LDA, E, IPIV, W, LDW, */
/*                             INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KB, LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), E( * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > SLASYF_RK computes a partial factorization of a real symmetric */
/* > matrix A using the bounded Bunch-Kaufman (rook) diagonal */
/* > pivoting method. The partial factorization has the form: */
/* > */
/* > A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or: */
/* >       ( 0  U22 ) (  0   D  ) ( U12**T U22**T ) */
/* > */
/* > A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L', */
/* >       ( L21  I ) (  0  A22 ) (  0       I    ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > */
/* > SLASYF_RK is an auxiliary routine called by SSYTRF_RK. It uses */
/* > blocked code (calling Level 3 BLAS) to update the submatrix */
/* > A11 (if UPLO = 'U') or A22 (if UPLO = 'L'). */
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
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The maximum number of columns of the matrix A that should be */
/* >          factored.  NB should be at least 2 to allow for 2-by-2 pivot */
/* >          blocks. */
/* > \endverbatim */
/* > */
/* > \param[out] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of columns of A that were actually factored. */
/* >          KB is either NB-1 or NB, or N if N <= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          E is REAL array, dimension (N) */
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
/* >          at each factorization step. */
/* > */
/* >          If UPLO = 'U', */
/* >          ( in factorization order, k decreases from N to 1 ): */
/* >            a) A single positive entry IPIV(k) > 0 means: */
/* >               D(k,k) is a 1-by-1 diagonal block. */
/* >               If IPIV(k) != k, rows and columns k and IPIV(k) were */
/* >               interchanged in the submatrix A(1:N,N-KB+1:N); */
/* >               If IPIV(k) = k, no interchange occurred. */
/* > */
/* > */
/* >            b) A pair of consecutive negative entries */
/* >               IPIV(k) < 0 and IPIV(k-1) < 0 means: */
/* >               D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
/* >               (NOTE: negative entries in IPIV appear ONLY in pairs). */
/* >               1) If -IPIV(k) != k, rows and columns */
/* >                  k and -IPIV(k) were interchanged */
/* >                  in the matrix A(1:N,N-KB+1:N). */
/* >                  If -IPIV(k) = k, no interchange occurred. */
/* >               2) If -IPIV(k-1) != k-1, rows and columns */
/* >                  k-1 and -IPIV(k-1) were interchanged */
/* >                  in the submatrix A(1:N,N-KB+1:N). */
/* >                  If -IPIV(k-1) = k-1, no interchange occurred. */
/* > */
/* >            c) In both cases a) and b) is always ABS( IPIV(k) ) <= k. */
/* > */
/* >            d) NOTE: Any entry IPIV(k) is always NONZERO on output. */
/* > */
/* >          If UPLO = 'L', */
/* >          ( in factorization order, k increases from 1 to N ): */
/* >            a) A single positive entry IPIV(k) > 0 means: */
/* >               D(k,k) is a 1-by-1 diagonal block. */
/* >               If IPIV(k) != k, rows and columns k and IPIV(k) were */
/* >               interchanged in the submatrix A(1:N,1:KB). */
/* >               If IPIV(k) = k, no interchange occurred. */
/* > */
/* >            b) A pair of consecutive negative entries */
/* >               IPIV(k) < 0 and IPIV(k+1) < 0 means: */
/* >               D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* >               (NOTE: negative entries in IPIV appear ONLY in pairs). */
/* >               1) If -IPIV(k) != k, rows and columns */
/* >                  k and -IPIV(k) were interchanged */
/* >                  in the submatrix A(1:N,1:KB). */
/* >                  If -IPIV(k) = k, no interchange occurred. */
/* >               2) If -IPIV(k+1) != k+1, rows and columns */
/* >                  k-1 and -IPIV(k-1) were interchanged */
/* >                  in the submatrix A(1:N,1:KB). */
/* >                  If -IPIV(k+1) = k+1, no interchange occurred. */
/* > */
/* >            c) In both cases a) and b) is always ABS( IPIV(k) ) >= k. */
/* > */
/* >            d) NOTE: Any entry IPIV(k) is always NONZERO on output. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (LDW,NB) */
/* > \endverbatim */
/* > */
/* > \param[in] LDW */
/* > \verbatim */
/* >          LDW is INTEGER */
/* >          The leading dimension of the array W.  LDW >= max(1,N). */
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

/* > \ingroup singleSYcomputational */

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
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int slasyf_rk__(char *uplo, integer *n, integer *nb, integer 
	*kb, doublereal *a, integer *lda, doublereal *e, integer *ipiv, 
	doublereal *w, integer *ldw, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, p;
    static doublereal t, r1, d11, d12, d21, d22;
    static integer jb, ii, jj, kk, kp, kw, kkw;
    static logical done;
    static integer imax, jmax;
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal sfmin;
    static integer itemp;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer kstep;
    static doublereal stemp;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal absakk;
    extern doublereal slamch_(char *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
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

#line 308 "slasyf_rk.f"
    /* Parameter adjustments */
#line 308 "slasyf_rk.f"
    a_dim1 = *lda;
#line 308 "slasyf_rk.f"
    a_offset = 1 + a_dim1;
#line 308 "slasyf_rk.f"
    a -= a_offset;
#line 308 "slasyf_rk.f"
    --e;
#line 308 "slasyf_rk.f"
    --ipiv;
#line 308 "slasyf_rk.f"
    w_dim1 = *ldw;
#line 308 "slasyf_rk.f"
    w_offset = 1 + w_dim1;
#line 308 "slasyf_rk.f"
    w -= w_offset;
#line 308 "slasyf_rk.f"

#line 308 "slasyf_rk.f"
    /* Function Body */
#line 308 "slasyf_rk.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 312 "slasyf_rk.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 316 "slasyf_rk.f"
    sfmin = slamch_("S", (ftnlen)1);

#line 318 "slasyf_rk.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 */

/*        Initilize the first entry of array E, where superdiagonal */
/*        elements of D are stored */

#line 327 "slasyf_rk.f"
	e[1] = 0.;

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

#line 331 "slasyf_rk.f"
	k = *n;
#line 332 "slasyf_rk.f"
L10:

/*        KW is the column of W which corresponds to column K of A */

#line 336 "slasyf_rk.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 340 "slasyf_rk.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 340 "slasyf_rk.f"
	    goto L30;
#line 340 "slasyf_rk.f"
	}

#line 343 "slasyf_rk.f"
	kstep = 1;
#line 344 "slasyf_rk.f"
	p = k;

/*        Copy column K of A to column KW of W and update it */

#line 348 "slasyf_rk.f"
	scopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
#line 349 "slasyf_rk.f"
	if (k < *n) {
#line 349 "slasyf_rk.f"
	    i__1 = *n - k;
#line 349 "slasyf_rk.f"
	    sgemv_("No transpose", &k, &i__1, &c_b9, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b10, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 349 "slasyf_rk.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 356 "slasyf_rk.f"
	absakk = (d__1 = w[k + kw * w_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 362 "slasyf_rk.f"
	if (k > 1) {
#line 363 "slasyf_rk.f"
	    i__1 = k - 1;
#line 363 "slasyf_rk.f"
	    imax = isamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 364 "slasyf_rk.f"
	    colmax = (d__1 = w[imax + kw * w_dim1], abs(d__1));
#line 365 "slasyf_rk.f"
	} else {
#line 366 "slasyf_rk.f"
	    colmax = 0.;
#line 367 "slasyf_rk.f"
	}

#line 369 "slasyf_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 373 "slasyf_rk.f"
	    if (*info == 0) {
#line 373 "slasyf_rk.f"
		*info = k;
#line 373 "slasyf_rk.f"
	    }
#line 375 "slasyf_rk.f"
	    kp = k;
#line 376 "slasyf_rk.f"
	    scopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);

/*           Set E( K ) to zero */

#line 380 "slasyf_rk.f"
	    if (k > 1) {
#line 380 "slasyf_rk.f"
		e[k] = 0.;
#line 380 "slasyf_rk.f"
	    }

#line 383 "slasyf_rk.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 392 "slasyf_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 396 "slasyf_rk.f"
		kp = k;

#line 398 "slasyf_rk.f"
	    } else {

#line 400 "slasyf_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 404 "slasyf_rk.f"
L12:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column KW-1 of W and update it */

#line 411 "slasyf_rk.f"
		scopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
#line 412 "slasyf_rk.f"
		i__1 = k - imax;
#line 412 "slasyf_rk.f"
		scopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);

#line 415 "slasyf_rk.f"
		if (k < *n) {
#line 415 "slasyf_rk.f"
		    i__1 = *n - k;
#line 415 "slasyf_rk.f"
		    sgemv_("No transpose", &k, &i__1, &c_b9, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b10, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 415 "slasyf_rk.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 424 "slasyf_rk.f"
		if (imax != k) {
#line 425 "slasyf_rk.f"
		    i__1 = k - imax;
#line 425 "slasyf_rk.f"
		    jmax = imax + isamax_(&i__1, &w[imax + 1 + (kw - 1) * 
			    w_dim1], &c__1);
#line 427 "slasyf_rk.f"
		    rowmax = (d__1 = w[jmax + (kw - 1) * w_dim1], abs(d__1));
#line 428 "slasyf_rk.f"
		} else {
#line 429 "slasyf_rk.f"
		    rowmax = 0.;
#line 430 "slasyf_rk.f"
		}

#line 432 "slasyf_rk.f"
		if (imax > 1) {
#line 433 "slasyf_rk.f"
		    i__1 = imax - 1;
#line 433 "slasyf_rk.f"
		    itemp = isamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
#line 434 "slasyf_rk.f"
		    stemp = (d__1 = w[itemp + (kw - 1) * w_dim1], abs(d__1));
#line 435 "slasyf_rk.f"
		    if (stemp > rowmax) {
#line 436 "slasyf_rk.f"
			rowmax = stemp;
#line 437 "slasyf_rk.f"
			jmax = itemp;
#line 438 "slasyf_rk.f"
		    }
#line 439 "slasyf_rk.f"
		}

/*                 Equivalent to testing for */
/*                 ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 445 "slasyf_rk.f"
		if (! ((d__1 = w[imax + (kw - 1) * w_dim1], abs(d__1)) < 
			alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 451 "slasyf_rk.f"
		    kp = imax;

/*                    copy column KW-1 of W to column KW of W */

#line 455 "slasyf_rk.f"
		    scopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 457 "slasyf_rk.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 462 "slasyf_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 468 "slasyf_rk.f"
		    kp = imax;
#line 469 "slasyf_rk.f"
		    kstep = 2;
#line 470 "slasyf_rk.f"
		    done = TRUE_;
#line 471 "slasyf_rk.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 475 "slasyf_rk.f"
		    p = imax;
#line 476 "slasyf_rk.f"
		    colmax = rowmax;
#line 477 "slasyf_rk.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 481 "slasyf_rk.f"
		    scopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 483 "slasyf_rk.f"
		}

/*                 End pivot search loop body */

#line 487 "slasyf_rk.f"
		if (! done) {
#line 487 "slasyf_rk.f"
		    goto L12;
#line 487 "slasyf_rk.f"
		}

#line 489 "slasyf_rk.f"
	    }

/*           ============================================================ */

#line 493 "slasyf_rk.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 497 "slasyf_rk.f"
	    kkw = *nb + kk - *n;

#line 499 "slasyf_rk.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

#line 503 "slasyf_rk.f"
		i__1 = k - p;
#line 503 "slasyf_rk.f"
		scopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * 
			a_dim1], lda);
#line 504 "slasyf_rk.f"
		scopy_(&p, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &
			c__1);

/*              Interchange rows K and P in last N-K+1 columns of A */
/*              and last N-K+2 columns of W */

#line 509 "slasyf_rk.f"
		i__1 = *n - k + 1;
#line 509 "slasyf_rk.f"
		sswap_(&i__1, &a[k + k * a_dim1], lda, &a[p + k * a_dim1], 
			lda);
#line 510 "slasyf_rk.f"
		i__1 = *n - kk + 1;
#line 510 "slasyf_rk.f"
		sswap_(&i__1, &w[k + kkw * w_dim1], ldw, &w[p + kkw * w_dim1],
			 ldw);
#line 511 "slasyf_rk.f"
	    }

/*           Updated column KP is already stored in column KKW of W */

#line 515 "slasyf_rk.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

#line 519 "slasyf_rk.f"
		a[kp + k * a_dim1] = a[kk + k * a_dim1];
#line 520 "slasyf_rk.f"
		i__1 = k - 1 - kp;
#line 520 "slasyf_rk.f"
		scopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 522 "slasyf_rk.f"
		scopy_(&kp, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
			c__1);

/*              Interchange rows KK and KP in last N-KK+1 columns */
/*              of A and W */

#line 527 "slasyf_rk.f"
		i__1 = *n - kk + 1;
#line 527 "slasyf_rk.f"
		sswap_(&i__1, &a[kk + kk * a_dim1], lda, &a[kp + kk * a_dim1],
			 lda);
#line 528 "slasyf_rk.f"
		i__1 = *n - kk + 1;
#line 528 "slasyf_rk.f"
		sswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 530 "slasyf_rk.f"
	    }

#line 532 "slasyf_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column KW of W now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Store U(k) in column k of A */

#line 542 "slasyf_rk.f"
		scopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 543 "slasyf_rk.f"
		if (k > 1) {
#line 544 "slasyf_rk.f"
		    if ((d__1 = a[k + k * a_dim1], abs(d__1)) >= sfmin) {
#line 545 "slasyf_rk.f"
			r1 = 1. / a[k + k * a_dim1];
#line 546 "slasyf_rk.f"
			i__1 = k - 1;
#line 546 "slasyf_rk.f"
			sscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 547 "slasyf_rk.f"
		    } else if (a[k + k * a_dim1] != 0.) {
#line 548 "slasyf_rk.f"
			i__1 = k - 1;
#line 548 "slasyf_rk.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 549 "slasyf_rk.f"
			    a[ii + k * a_dim1] /= a[k + k * a_dim1];
#line 550 "slasyf_rk.f"
/* L14: */
#line 550 "slasyf_rk.f"
			}
#line 551 "slasyf_rk.f"
		    }

/*                 Store the superdiagonal element of D in array E */

#line 555 "slasyf_rk.f"
		    e[k] = 0.;

#line 557 "slasyf_rk.f"
		}

#line 559 "slasyf_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns KW and KW-1 of W now */
/*              hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

#line 569 "slasyf_rk.f"
		if (k > 2) {

/*                 Store U(k) and U(k-1) in columns k and k-1 of A */

#line 573 "slasyf_rk.f"
		    d12 = w[k - 1 + kw * w_dim1];
#line 574 "slasyf_rk.f"
		    d11 = w[k + kw * w_dim1] / d12;
#line 575 "slasyf_rk.f"
		    d22 = w[k - 1 + (kw - 1) * w_dim1] / d12;
#line 576 "slasyf_rk.f"
		    t = 1. / (d11 * d22 - 1.);
#line 577 "slasyf_rk.f"
		    i__1 = k - 2;
#line 577 "slasyf_rk.f"
		    for (j = 1; j <= i__1; ++j) {
#line 578 "slasyf_rk.f"
			a[j + (k - 1) * a_dim1] = t * ((d11 * w[j + (kw - 1) *
				 w_dim1] - w[j + kw * w_dim1]) / d12);
#line 580 "slasyf_rk.f"
			a[j + k * a_dim1] = t * ((d22 * w[j + kw * w_dim1] - 
				w[j + (kw - 1) * w_dim1]) / d12);
#line 582 "slasyf_rk.f"
/* L20: */
#line 582 "slasyf_rk.f"
		    }
#line 583 "slasyf_rk.f"
		}

/*              Copy diagonal elements of D(K) to A, */
/*              copy superdiagonal element of D(K) to E(K) and */
/*              ZERO out superdiagonal entry of A */

#line 589 "slasyf_rk.f"
		a[k - 1 + (k - 1) * a_dim1] = w[k - 1 + (kw - 1) * w_dim1];
#line 590 "slasyf_rk.f"
		a[k - 1 + k * a_dim1] = 0.;
#line 591 "slasyf_rk.f"
		a[k + k * a_dim1] = w[k + kw * w_dim1];
#line 592 "slasyf_rk.f"
		e[k] = w[k - 1 + kw * w_dim1];
#line 593 "slasyf_rk.f"
		e[k - 1] = 0.;

#line 595 "slasyf_rk.f"
	    }

/*           End column K is nonsingular */

#line 599 "slasyf_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 603 "slasyf_rk.f"
	if (kstep == 1) {
#line 604 "slasyf_rk.f"
	    ipiv[k] = kp;
#line 605 "slasyf_rk.f"
	} else {
#line 606 "slasyf_rk.f"
	    ipiv[k] = -p;
#line 607 "slasyf_rk.f"
	    ipiv[k - 1] = -kp;
#line 608 "slasyf_rk.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 612 "slasyf_rk.f"
	k -= kstep;
#line 613 "slasyf_rk.f"
	goto L10;

#line 615 "slasyf_rk.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T */

/*        computing blocks of NB columns at a time */

#line 623 "slasyf_rk.f"
	i__1 = -(*nb);
#line 623 "slasyf_rk.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 624 "slasyf_rk.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 624 "slasyf_rk.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 628 "slasyf_rk.f"
	    i__2 = j + jb - 1;
#line 628 "slasyf_rk.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 629 "slasyf_rk.f"
		i__3 = jj - j + 1;
#line 629 "slasyf_rk.f"
		i__4 = *n - k;
#line 629 "slasyf_rk.f"
		sgemv_("No transpose", &i__3, &i__4, &c_b9, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b10,
			 &a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 632 "slasyf_rk.f"
/* L40: */
#line 632 "slasyf_rk.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 636 "slasyf_rk.f"
	    if (j >= 2) {
#line 636 "slasyf_rk.f"
		i__2 = j - 1;
#line 636 "slasyf_rk.f"
		i__3 = *n - k;
#line 636 "slasyf_rk.f"
		sgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &c_b9, 
			&a[(k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * 
			w_dim1], ldw, &c_b10, &a[j * a_dim1 + 1], lda, (
			ftnlen)12, (ftnlen)9);
#line 636 "slasyf_rk.f"
	    }
#line 640 "slasyf_rk.f"
/* L50: */
#line 640 "slasyf_rk.f"
	}

/*        Set KB to the number of columns factorized */

#line 644 "slasyf_rk.f"
	*kb = *n - k;

#line 646 "slasyf_rk.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 */

/*        Initilize the unused last entry of the subdiagonal array E. */

#line 654 "slasyf_rk.f"
	e[*n] = 0.;

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 658 "slasyf_rk.f"
	k = 1;
#line 659 "slasyf_rk.f"
L70:

/*        Exit from loop */

#line 663 "slasyf_rk.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 663 "slasyf_rk.f"
	    goto L90;
#line 663 "slasyf_rk.f"
	}

#line 666 "slasyf_rk.f"
	kstep = 1;
#line 667 "slasyf_rk.f"
	p = k;

/*        Copy column K of A to column K of W and update it */

#line 671 "slasyf_rk.f"
	i__1 = *n - k + 1;
#line 671 "slasyf_rk.f"
	scopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
#line 672 "slasyf_rk.f"
	if (k > 1) {
#line 672 "slasyf_rk.f"
	    i__1 = *n - k + 1;
#line 672 "slasyf_rk.f"
	    i__2 = k - 1;
#line 672 "slasyf_rk.f"
	    sgemv_("No transpose", &i__1, &i__2, &c_b9, &a[k + a_dim1], lda, &
		    w[k + w_dim1], ldw, &c_b10, &w[k + k * w_dim1], &c__1, (
		    ftnlen)12);
#line 672 "slasyf_rk.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 679 "slasyf_rk.f"
	absakk = (d__1 = w[k + k * w_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 685 "slasyf_rk.f"
	if (k < *n) {
#line 686 "slasyf_rk.f"
	    i__1 = *n - k;
#line 686 "slasyf_rk.f"
	    imax = k + isamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 687 "slasyf_rk.f"
	    colmax = (d__1 = w[imax + k * w_dim1], abs(d__1));
#line 688 "slasyf_rk.f"
	} else {
#line 689 "slasyf_rk.f"
	    colmax = 0.;
#line 690 "slasyf_rk.f"
	}

#line 692 "slasyf_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 696 "slasyf_rk.f"
	    if (*info == 0) {
#line 696 "slasyf_rk.f"
		*info = k;
#line 696 "slasyf_rk.f"
	    }
#line 698 "slasyf_rk.f"
	    kp = k;
#line 699 "slasyf_rk.f"
	    i__1 = *n - k + 1;
#line 699 "slasyf_rk.f"
	    scopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
		    c__1);

/*           Set E( K ) to zero */

#line 703 "slasyf_rk.f"
	    if (k < *n) {
#line 703 "slasyf_rk.f"
		e[k] = 0.;
#line 703 "slasyf_rk.f"
	    }

#line 706 "slasyf_rk.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 715 "slasyf_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 719 "slasyf_rk.f"
		kp = k;

#line 721 "slasyf_rk.f"
	    } else {

#line 723 "slasyf_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 727 "slasyf_rk.f"
L72:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column K+1 of W and update it */

#line 734 "slasyf_rk.f"
		i__1 = imax - k;
#line 734 "slasyf_rk.f"
		scopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 735 "slasyf_rk.f"
		i__1 = *n - imax + 1;
#line 735 "slasyf_rk.f"
		scopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 
			1) * w_dim1], &c__1);
#line 737 "slasyf_rk.f"
		if (k > 1) {
#line 737 "slasyf_rk.f"
		    i__1 = *n - k + 1;
#line 737 "slasyf_rk.f"
		    i__2 = k - 1;
#line 737 "slasyf_rk.f"
		    sgemv_("No transpose", &i__1, &i__2, &c_b9, &a[k + a_dim1]
			    , lda, &w[imax + w_dim1], ldw, &c_b10, &w[k + (k 
			    + 1) * w_dim1], &c__1, (ftnlen)12);
#line 737 "slasyf_rk.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 746 "slasyf_rk.f"
		if (imax != k) {
#line 747 "slasyf_rk.f"
		    i__1 = imax - k;
#line 747 "slasyf_rk.f"
		    jmax = k - 1 + isamax_(&i__1, &w[k + (k + 1) * w_dim1], &
			    c__1);
#line 748 "slasyf_rk.f"
		    rowmax = (d__1 = w[jmax + (k + 1) * w_dim1], abs(d__1));
#line 749 "slasyf_rk.f"
		} else {
#line 750 "slasyf_rk.f"
		    rowmax = 0.;
#line 751 "slasyf_rk.f"
		}

#line 753 "slasyf_rk.f"
		if (imax < *n) {
#line 754 "slasyf_rk.f"
		    i__1 = *n - imax;
#line 754 "slasyf_rk.f"
		    itemp = imax + isamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
#line 755 "slasyf_rk.f"
		    stemp = (d__1 = w[itemp + (k + 1) * w_dim1], abs(d__1));
#line 756 "slasyf_rk.f"
		    if (stemp > rowmax) {
#line 757 "slasyf_rk.f"
			rowmax = stemp;
#line 758 "slasyf_rk.f"
			jmax = itemp;
#line 759 "slasyf_rk.f"
		    }
#line 760 "slasyf_rk.f"
		}

/*                 Equivalent to testing for */
/*                 ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 766 "slasyf_rk.f"
		if (! ((d__1 = w[imax + (k + 1) * w_dim1], abs(d__1)) < alpha 
			* rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 772 "slasyf_rk.f"
		    kp = imax;

/*                    copy column K+1 of W to column K of W */

#line 776 "slasyf_rk.f"
		    i__1 = *n - k + 1;
#line 776 "slasyf_rk.f"
		    scopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 778 "slasyf_rk.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 783 "slasyf_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 789 "slasyf_rk.f"
		    kp = imax;
#line 790 "slasyf_rk.f"
		    kstep = 2;
#line 791 "slasyf_rk.f"
		    done = TRUE_;
#line 792 "slasyf_rk.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 796 "slasyf_rk.f"
		    p = imax;
#line 797 "slasyf_rk.f"
		    colmax = rowmax;
#line 798 "slasyf_rk.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 802 "slasyf_rk.f"
		    i__1 = *n - k + 1;
#line 802 "slasyf_rk.f"
		    scopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 804 "slasyf_rk.f"
		}

/*                 End pivot search loop body */

#line 808 "slasyf_rk.f"
		if (! done) {
#line 808 "slasyf_rk.f"
		    goto L72;
#line 808 "slasyf_rk.f"
		}

#line 810 "slasyf_rk.f"
	    }

/*           ============================================================ */

#line 814 "slasyf_rk.f"
	    kk = k + kstep - 1;

#line 816 "slasyf_rk.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

#line 820 "slasyf_rk.f"
		i__1 = p - k;
#line 820 "slasyf_rk.f"
		scopy_(&i__1, &a[k + k * a_dim1], &c__1, &a[p + k * a_dim1], 
			lda);
#line 821 "slasyf_rk.f"
		i__1 = *n - p + 1;
#line 821 "slasyf_rk.f"
		scopy_(&i__1, &a[p + k * a_dim1], &c__1, &a[p + p * a_dim1], &
			c__1);

/*              Interchange rows K and P in first K columns of A */
/*              and first K+1 columns of W */

#line 826 "slasyf_rk.f"
		sswap_(&k, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 827 "slasyf_rk.f"
		sswap_(&kk, &w[k + w_dim1], ldw, &w[p + w_dim1], ldw);
#line 828 "slasyf_rk.f"
	    }

/*           Updated column KP is already stored in column KK of W */

#line 832 "slasyf_rk.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

#line 836 "slasyf_rk.f"
		a[kp + k * a_dim1] = a[kk + k * a_dim1];
#line 837 "slasyf_rk.f"
		i__1 = kp - k - 1;
#line 837 "slasyf_rk.f"
		scopy_(&i__1, &a[k + 1 + kk * a_dim1], &c__1, &a[kp + (k + 1) 
			* a_dim1], lda);
#line 838 "slasyf_rk.f"
		i__1 = *n - kp + 1;
#line 838 "slasyf_rk.f"
		scopy_(&i__1, &a[kp + kk * a_dim1], &c__1, &a[kp + kp * 
			a_dim1], &c__1);

/*              Interchange rows KK and KP in first KK columns of A and W */

#line 842 "slasyf_rk.f"
		sswap_(&kk, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 843 "slasyf_rk.f"
		sswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 844 "slasyf_rk.f"
	    }

#line 846 "slasyf_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of W now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

/*              Store L(k) in column k of A */

#line 856 "slasyf_rk.f"
		i__1 = *n - k + 1;
#line 856 "slasyf_rk.f"
		scopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 857 "slasyf_rk.f"
		if (k < *n) {
#line 858 "slasyf_rk.f"
		    if ((d__1 = a[k + k * a_dim1], abs(d__1)) >= sfmin) {
#line 859 "slasyf_rk.f"
			r1 = 1. / a[k + k * a_dim1];
#line 860 "slasyf_rk.f"
			i__1 = *n - k;
#line 860 "slasyf_rk.f"
			sscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 861 "slasyf_rk.f"
		    } else if (a[k + k * a_dim1] != 0.) {
#line 862 "slasyf_rk.f"
			i__1 = *n;
#line 862 "slasyf_rk.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 863 "slasyf_rk.f"
			    a[ii + k * a_dim1] /= a[k + k * a_dim1];
#line 864 "slasyf_rk.f"
/* L74: */
#line 864 "slasyf_rk.f"
			}
#line 865 "slasyf_rk.f"
		    }

/*                 Store the subdiagonal element of D in array E */

#line 869 "slasyf_rk.f"
		    e[k] = 0.;

#line 871 "slasyf_rk.f"
		}

#line 873 "slasyf_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

#line 882 "slasyf_rk.f"
		if (k < *n - 1) {

/*                 Store L(k) and L(k+1) in columns k and k+1 of A */

#line 886 "slasyf_rk.f"
		    d21 = w[k + 1 + k * w_dim1];
#line 887 "slasyf_rk.f"
		    d11 = w[k + 1 + (k + 1) * w_dim1] / d21;
#line 888 "slasyf_rk.f"
		    d22 = w[k + k * w_dim1] / d21;
#line 889 "slasyf_rk.f"
		    t = 1. / (d11 * d22 - 1.);
#line 890 "slasyf_rk.f"
		    i__1 = *n;
#line 890 "slasyf_rk.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 891 "slasyf_rk.f"
			a[j + k * a_dim1] = t * ((d11 * w[j + k * w_dim1] - w[
				j + (k + 1) * w_dim1]) / d21);
#line 893 "slasyf_rk.f"
			a[j + (k + 1) * a_dim1] = t * ((d22 * w[j + (k + 1) * 
				w_dim1] - w[j + k * w_dim1]) / d21);
#line 895 "slasyf_rk.f"
/* L80: */
#line 895 "slasyf_rk.f"
		    }
#line 896 "slasyf_rk.f"
		}

/*              Copy diagonal elements of D(K) to A, */
/*              copy subdiagonal element of D(K) to E(K) and */
/*              ZERO out subdiagonal entry of A */

#line 902 "slasyf_rk.f"
		a[k + k * a_dim1] = w[k + k * w_dim1];
#line 903 "slasyf_rk.f"
		a[k + 1 + k * a_dim1] = 0.;
#line 904 "slasyf_rk.f"
		a[k + 1 + (k + 1) * a_dim1] = w[k + 1 + (k + 1) * w_dim1];
#line 905 "slasyf_rk.f"
		e[k] = w[k + 1 + k * w_dim1];
#line 906 "slasyf_rk.f"
		e[k + 1] = 0.;

#line 908 "slasyf_rk.f"
	    }

/*           End column K is nonsingular */

#line 912 "slasyf_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 916 "slasyf_rk.f"
	if (kstep == 1) {
#line 917 "slasyf_rk.f"
	    ipiv[k] = kp;
#line 918 "slasyf_rk.f"
	} else {
#line 919 "slasyf_rk.f"
	    ipiv[k] = -p;
#line 920 "slasyf_rk.f"
	    ipiv[k + 1] = -kp;
#line 921 "slasyf_rk.f"
	}

/*        Increase K and return to the start of the main loop */

#line 925 "slasyf_rk.f"
	k += kstep;
#line 926 "slasyf_rk.f"
	goto L70;

#line 928 "slasyf_rk.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T */

/*        computing blocks of NB columns at a time */

#line 936 "slasyf_rk.f"
	i__1 = *n;
#line 936 "slasyf_rk.f"
	i__2 = *nb;
#line 936 "slasyf_rk.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 937 "slasyf_rk.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 937 "slasyf_rk.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 941 "slasyf_rk.f"
	    i__3 = j + jb - 1;
#line 941 "slasyf_rk.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 942 "slasyf_rk.f"
		i__4 = j + jb - jj;
#line 942 "slasyf_rk.f"
		i__5 = k - 1;
#line 942 "slasyf_rk.f"
		sgemv_("No transpose", &i__4, &i__5, &c_b9, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b10, &a[jj + jj * 
			a_dim1], &c__1, (ftnlen)12);
#line 945 "slasyf_rk.f"
/* L100: */
#line 945 "slasyf_rk.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 949 "slasyf_rk.f"
	    if (j + jb <= *n) {
#line 949 "slasyf_rk.f"
		i__3 = *n - j - jb + 1;
#line 949 "slasyf_rk.f"
		i__4 = k - 1;
#line 949 "slasyf_rk.f"
		sgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &c_b9, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b10,
			 &a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 949 "slasyf_rk.f"
	    }
#line 953 "slasyf_rk.f"
/* L110: */
#line 953 "slasyf_rk.f"
	}

/*        Set KB to the number of columns factorized */

#line 957 "slasyf_rk.f"
	*kb = k - 1;

#line 959 "slasyf_rk.f"
    }

#line 961 "slasyf_rk.f"
    return 0;

/*     End of SLASYF_RK */

} /* slasyf_rk__ */

