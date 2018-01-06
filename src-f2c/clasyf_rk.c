#line 1 "clasyf_rk.f"
/* clasyf_rk.f -- translated by f2c (version 20100827).
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

#line 1 "clasyf_rk.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLASYF_RK computes a partial factorization of a complex symmetric indefinite matrix using bound
ed Bunch-Kaufman (rook) diagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLASYF_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasyf_
rk.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasyf_
rk.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasyf_
rk.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLASYF_RK( UPLO, N, NB, KB, A, LDA, E, IPIV, W, LDW, */
/*                             INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KB, LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), E( * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > CLASYF_RK computes a partial factorization of a complex symmetric */
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
/* > CLASYF_RK is an auxiliary routine called by CSYTRF_RK. It uses */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          E is COMPLEX array, dimension (N) */
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
/* >          W is COMPLEX array, dimension (LDW,NB) */
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

/* > \ingroup complexSYcomputational */

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
/* Subroutine */ int clasyf_rk__(char *uplo, integer *n, integer *nb, integer 
	*kb, doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	doublecomplex *w, integer *ldw, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, p;
    static doublecomplex t, r1, d11, d12, d21, d22;
    static integer jb, ii, jj, kk, kp, kw, kkw;
    static logical done;
    static integer imax, jmax;
    static doublereal alpha;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal sfmin;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer itemp;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer kstep;
    static doublereal stemp, absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 317 "clasyf_rk.f"
    /* Parameter adjustments */
#line 317 "clasyf_rk.f"
    a_dim1 = *lda;
#line 317 "clasyf_rk.f"
    a_offset = 1 + a_dim1;
#line 317 "clasyf_rk.f"
    a -= a_offset;
#line 317 "clasyf_rk.f"
    --e;
#line 317 "clasyf_rk.f"
    --ipiv;
#line 317 "clasyf_rk.f"
    w_dim1 = *ldw;
#line 317 "clasyf_rk.f"
    w_offset = 1 + w_dim1;
#line 317 "clasyf_rk.f"
    w -= w_offset;
#line 317 "clasyf_rk.f"

#line 317 "clasyf_rk.f"
    /* Function Body */
#line 317 "clasyf_rk.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 321 "clasyf_rk.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 325 "clasyf_rk.f"
    sfmin = slamch_("S", (ftnlen)1);

#line 327 "clasyf_rk.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 */

/*        Initilize the first entry of array E, where superdiagonal */
/*        elements of D are stored */

#line 336 "clasyf_rk.f"
	e[1].r = 0., e[1].i = 0.;

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

#line 340 "clasyf_rk.f"
	k = *n;
#line 341 "clasyf_rk.f"
L10:

/*        KW is the column of W which corresponds to column K of A */

#line 345 "clasyf_rk.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 349 "clasyf_rk.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 349 "clasyf_rk.f"
	    goto L30;
#line 349 "clasyf_rk.f"
	}

#line 352 "clasyf_rk.f"
	kstep = 1;
#line 353 "clasyf_rk.f"
	p = k;

/*        Copy column K of A to column KW of W and update it */

#line 357 "clasyf_rk.f"
	ccopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
#line 358 "clasyf_rk.f"
	if (k < *n) {
#line 358 "clasyf_rk.f"
	    i__1 = *n - k;
#line 358 "clasyf_rk.f"
	    z__1.r = -1., z__1.i = -0.;
#line 358 "clasyf_rk.f"
	    cgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 358 "clasyf_rk.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 365 "clasyf_rk.f"
	i__1 = k + kw * w_dim1;
#line 365 "clasyf_rk.f"
	absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[k + kw * 
		w_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 371 "clasyf_rk.f"
	if (k > 1) {
#line 372 "clasyf_rk.f"
	    i__1 = k - 1;
#line 372 "clasyf_rk.f"
	    imax = icamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 373 "clasyf_rk.f"
	    i__1 = imax + kw * w_dim1;
#line 373 "clasyf_rk.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    kw * w_dim1]), abs(d__2));
#line 374 "clasyf_rk.f"
	} else {
#line 375 "clasyf_rk.f"
	    colmax = 0.;
#line 376 "clasyf_rk.f"
	}

#line 378 "clasyf_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 382 "clasyf_rk.f"
	    if (*info == 0) {
#line 382 "clasyf_rk.f"
		*info = k;
#line 382 "clasyf_rk.f"
	    }
#line 384 "clasyf_rk.f"
	    kp = k;
#line 385 "clasyf_rk.f"
	    ccopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);

/*           Set E( K ) to zero */

#line 389 "clasyf_rk.f"
	    if (k > 1) {
#line 389 "clasyf_rk.f"
		i__1 = k;
#line 389 "clasyf_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 389 "clasyf_rk.f"
	    }

#line 392 "clasyf_rk.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 401 "clasyf_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 405 "clasyf_rk.f"
		kp = k;

#line 407 "clasyf_rk.f"
	    } else {

#line 409 "clasyf_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 413 "clasyf_rk.f"
L12:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column KW-1 of W and update it */

#line 420 "clasyf_rk.f"
		ccopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
#line 421 "clasyf_rk.f"
		i__1 = k - imax;
#line 421 "clasyf_rk.f"
		ccopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);

#line 424 "clasyf_rk.f"
		if (k < *n) {
#line 424 "clasyf_rk.f"
		    i__1 = *n - k;
#line 424 "clasyf_rk.f"
		    z__1.r = -1., z__1.i = -0.;
#line 424 "clasyf_rk.f"
		    cgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 424 "clasyf_rk.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 433 "clasyf_rk.f"
		if (imax != k) {
#line 434 "clasyf_rk.f"
		    i__1 = k - imax;
#line 434 "clasyf_rk.f"
		    jmax = imax + icamax_(&i__1, &w[imax + 1 + (kw - 1) * 
			    w_dim1], &c__1);
#line 436 "clasyf_rk.f"
		    i__1 = jmax + (kw - 1) * w_dim1;
#line 436 "clasyf_rk.f"
		    rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    w[jmax + (kw - 1) * w_dim1]), abs(d__2));
#line 437 "clasyf_rk.f"
		} else {
#line 438 "clasyf_rk.f"
		    rowmax = 0.;
#line 439 "clasyf_rk.f"
		}

#line 441 "clasyf_rk.f"
		if (imax > 1) {
#line 442 "clasyf_rk.f"
		    i__1 = imax - 1;
#line 442 "clasyf_rk.f"
		    itemp = icamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
#line 443 "clasyf_rk.f"
		    i__1 = itemp + (kw - 1) * w_dim1;
#line 443 "clasyf_rk.f"
		    stemp = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    itemp + (kw - 1) * w_dim1]), abs(d__2));
#line 444 "clasyf_rk.f"
		    if (stemp > rowmax) {
#line 445 "clasyf_rk.f"
			rowmax = stemp;
#line 446 "clasyf_rk.f"
			jmax = itemp;
#line 447 "clasyf_rk.f"
		    }
#line 448 "clasyf_rk.f"
		}

/*                 Equivalent to testing for */
/*                 CABS1( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 454 "clasyf_rk.f"
		i__1 = imax + (kw - 1) * w_dim1;
#line 454 "clasyf_rk.f"
		if (! ((d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax 
			+ (kw - 1) * w_dim1]), abs(d__2)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 460 "clasyf_rk.f"
		    kp = imax;

/*                    copy column KW-1 of W to column KW of W */

#line 464 "clasyf_rk.f"
		    ccopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 466 "clasyf_rk.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 471 "clasyf_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 477 "clasyf_rk.f"
		    kp = imax;
#line 478 "clasyf_rk.f"
		    kstep = 2;
#line 479 "clasyf_rk.f"
		    done = TRUE_;
#line 480 "clasyf_rk.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 484 "clasyf_rk.f"
		    p = imax;
#line 485 "clasyf_rk.f"
		    colmax = rowmax;
#line 486 "clasyf_rk.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 490 "clasyf_rk.f"
		    ccopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 492 "clasyf_rk.f"
		}

/*                 End pivot search loop body */

#line 496 "clasyf_rk.f"
		if (! done) {
#line 496 "clasyf_rk.f"
		    goto L12;
#line 496 "clasyf_rk.f"
		}

#line 498 "clasyf_rk.f"
	    }

/*           ============================================================ */

#line 502 "clasyf_rk.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 506 "clasyf_rk.f"
	    kkw = *nb + kk - *n;

#line 508 "clasyf_rk.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

#line 512 "clasyf_rk.f"
		i__1 = k - p;
#line 512 "clasyf_rk.f"
		ccopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * 
			a_dim1], lda);
#line 513 "clasyf_rk.f"
		ccopy_(&p, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &
			c__1);

/*              Interchange rows K and P in last N-K+1 columns of A */
/*              and last N-K+2 columns of W */

#line 518 "clasyf_rk.f"
		i__1 = *n - k + 1;
#line 518 "clasyf_rk.f"
		cswap_(&i__1, &a[k + k * a_dim1], lda, &a[p + k * a_dim1], 
			lda);
#line 519 "clasyf_rk.f"
		i__1 = *n - kk + 1;
#line 519 "clasyf_rk.f"
		cswap_(&i__1, &w[k + kkw * w_dim1], ldw, &w[p + kkw * w_dim1],
			 ldw);
#line 520 "clasyf_rk.f"
	    }

/*           Updated column KP is already stored in column KKW of W */

#line 524 "clasyf_rk.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

#line 528 "clasyf_rk.f"
		i__1 = kp + k * a_dim1;
#line 528 "clasyf_rk.f"
		i__2 = kk + k * a_dim1;
#line 528 "clasyf_rk.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 529 "clasyf_rk.f"
		i__1 = k - 1 - kp;
#line 529 "clasyf_rk.f"
		ccopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 531 "clasyf_rk.f"
		ccopy_(&kp, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
			c__1);

/*              Interchange rows KK and KP in last N-KK+1 columns */
/*              of A and W */

#line 536 "clasyf_rk.f"
		i__1 = *n - kk + 1;
#line 536 "clasyf_rk.f"
		cswap_(&i__1, &a[kk + kk * a_dim1], lda, &a[kp + kk * a_dim1],
			 lda);
#line 537 "clasyf_rk.f"
		i__1 = *n - kk + 1;
#line 537 "clasyf_rk.f"
		cswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 539 "clasyf_rk.f"
	    }

#line 541 "clasyf_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column KW of W now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Store U(k) in column k of A */

#line 551 "clasyf_rk.f"
		ccopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 552 "clasyf_rk.f"
		if (k > 1) {
#line 553 "clasyf_rk.f"
		    i__1 = k + k * a_dim1;
#line 553 "clasyf_rk.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + 
			    k * a_dim1]), abs(d__2)) >= sfmin) {
#line 554 "clasyf_rk.f"
			z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 554 "clasyf_rk.f"
			r1.r = z__1.r, r1.i = z__1.i;
#line 555 "clasyf_rk.f"
			i__1 = k - 1;
#line 555 "clasyf_rk.f"
			cscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 556 "clasyf_rk.f"
		    } else /* if(complicated condition) */ {
#line 556 "clasyf_rk.f"
			i__1 = k + k * a_dim1;
#line 556 "clasyf_rk.f"
			if (a[i__1].r != 0. || a[i__1].i != 0.) {
#line 557 "clasyf_rk.f"
			    i__1 = k - 1;
#line 557 "clasyf_rk.f"
			    for (ii = 1; ii <= i__1; ++ii) {
#line 558 "clasyf_rk.f"
				i__2 = ii + k * a_dim1;
#line 558 "clasyf_rk.f"
				z_div(&z__1, &a[ii + k * a_dim1], &a[k + k * 
					a_dim1]);
#line 558 "clasyf_rk.f"
				a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 559 "clasyf_rk.f"
/* L14: */
#line 559 "clasyf_rk.f"
			    }
#line 560 "clasyf_rk.f"
			}
#line 560 "clasyf_rk.f"
		    }

/*                 Store the superdiagonal element of D in array E */

#line 564 "clasyf_rk.f"
		    i__1 = k;
#line 564 "clasyf_rk.f"
		    e[i__1].r = 0., e[i__1].i = 0.;

#line 566 "clasyf_rk.f"
		}

#line 568 "clasyf_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns KW and KW-1 of W now */
/*              hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

#line 578 "clasyf_rk.f"
		if (k > 2) {

/*                 Store U(k) and U(k-1) in columns k and k-1 of A */

#line 582 "clasyf_rk.f"
		    i__1 = k - 1 + kw * w_dim1;
#line 582 "clasyf_rk.f"
		    d12.r = w[i__1].r, d12.i = w[i__1].i;
#line 583 "clasyf_rk.f"
		    z_div(&z__1, &w[k + kw * w_dim1], &d12);
#line 583 "clasyf_rk.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 584 "clasyf_rk.f"
		    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d12);
#line 584 "clasyf_rk.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 585 "clasyf_rk.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 585 "clasyf_rk.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 585 "clasyf_rk.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 585 "clasyf_rk.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 586 "clasyf_rk.f"
		    i__1 = k - 2;
#line 586 "clasyf_rk.f"
		    for (j = 1; j <= i__1; ++j) {
#line 587 "clasyf_rk.f"
			i__2 = j + (k - 1) * a_dim1;
#line 587 "clasyf_rk.f"
			i__3 = j + (kw - 1) * w_dim1;
#line 587 "clasyf_rk.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 587 "clasyf_rk.f"
			i__4 = j + kw * w_dim1;
#line 587 "clasyf_rk.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 587 "clasyf_rk.f"
			z_div(&z__2, &z__3, &d12);
#line 587 "clasyf_rk.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 587 "clasyf_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 589 "clasyf_rk.f"
			i__2 = j + k * a_dim1;
#line 589 "clasyf_rk.f"
			i__3 = j + kw * w_dim1;
#line 589 "clasyf_rk.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 589 "clasyf_rk.f"
			i__4 = j + (kw - 1) * w_dim1;
#line 589 "clasyf_rk.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 589 "clasyf_rk.f"
			z_div(&z__2, &z__3, &d12);
#line 589 "clasyf_rk.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 589 "clasyf_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 591 "clasyf_rk.f"
/* L20: */
#line 591 "clasyf_rk.f"
		    }
#line 592 "clasyf_rk.f"
		}

/*              Copy diagonal elements of D(K) to A, */
/*              copy superdiagonal element of D(K) to E(K) and */
/*              ZERO out superdiagonal entry of A */

#line 598 "clasyf_rk.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 598 "clasyf_rk.f"
		i__2 = k - 1 + (kw - 1) * w_dim1;
#line 598 "clasyf_rk.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 599 "clasyf_rk.f"
		i__1 = k - 1 + k * a_dim1;
#line 599 "clasyf_rk.f"
		a[i__1].r = 0., a[i__1].i = 0.;
#line 600 "clasyf_rk.f"
		i__1 = k + k * a_dim1;
#line 600 "clasyf_rk.f"
		i__2 = k + kw * w_dim1;
#line 600 "clasyf_rk.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 601 "clasyf_rk.f"
		i__1 = k;
#line 601 "clasyf_rk.f"
		i__2 = k - 1 + kw * w_dim1;
#line 601 "clasyf_rk.f"
		e[i__1].r = w[i__2].r, e[i__1].i = w[i__2].i;
#line 602 "clasyf_rk.f"
		i__1 = k - 1;
#line 602 "clasyf_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;

#line 604 "clasyf_rk.f"
	    }

/*           End column K is nonsingular */

#line 608 "clasyf_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 612 "clasyf_rk.f"
	if (kstep == 1) {
#line 613 "clasyf_rk.f"
	    ipiv[k] = kp;
#line 614 "clasyf_rk.f"
	} else {
#line 615 "clasyf_rk.f"
	    ipiv[k] = -p;
#line 616 "clasyf_rk.f"
	    ipiv[k - 1] = -kp;
#line 617 "clasyf_rk.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 621 "clasyf_rk.f"
	k -= kstep;
#line 622 "clasyf_rk.f"
	goto L10;

#line 624 "clasyf_rk.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T */

/*        computing blocks of NB columns at a time */

#line 632 "clasyf_rk.f"
	i__1 = -(*nb);
#line 632 "clasyf_rk.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 633 "clasyf_rk.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 633 "clasyf_rk.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 637 "clasyf_rk.f"
	    i__2 = j + jb - 1;
#line 637 "clasyf_rk.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 638 "clasyf_rk.f"
		i__3 = jj - j + 1;
#line 638 "clasyf_rk.f"
		i__4 = *n - k;
#line 638 "clasyf_rk.f"
		z__1.r = -1., z__1.i = -0.;
#line 638 "clasyf_rk.f"
		cgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, 
			&a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 641 "clasyf_rk.f"
/* L40: */
#line 641 "clasyf_rk.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 645 "clasyf_rk.f"
	    if (j >= 2) {
#line 645 "clasyf_rk.f"
		i__2 = j - 1;
#line 645 "clasyf_rk.f"
		i__3 = *n - k;
#line 645 "clasyf_rk.f"
		z__1.r = -1., z__1.i = -0.;
#line 645 "clasyf_rk.f"
		cgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, 
			&a[(k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * 
			w_dim1], ldw, &c_b1, &a[j * a_dim1 + 1], lda, (ftnlen)
			12, (ftnlen)9);
#line 645 "clasyf_rk.f"
	    }
#line 649 "clasyf_rk.f"
/* L50: */
#line 649 "clasyf_rk.f"
	}

/*        Set KB to the number of columns factorized */

#line 653 "clasyf_rk.f"
	*kb = *n - k;

#line 655 "clasyf_rk.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 */

/*        Initilize the unused last entry of the subdiagonal array E. */

#line 663 "clasyf_rk.f"
	i__1 = *n;
#line 663 "clasyf_rk.f"
	e[i__1].r = 0., e[i__1].i = 0.;

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 667 "clasyf_rk.f"
	k = 1;
#line 668 "clasyf_rk.f"
L70:

/*        Exit from loop */

#line 672 "clasyf_rk.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 672 "clasyf_rk.f"
	    goto L90;
#line 672 "clasyf_rk.f"
	}

#line 675 "clasyf_rk.f"
	kstep = 1;
#line 676 "clasyf_rk.f"
	p = k;

/*        Copy column K of A to column K of W and update it */

#line 680 "clasyf_rk.f"
	i__1 = *n - k + 1;
#line 680 "clasyf_rk.f"
	ccopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
#line 681 "clasyf_rk.f"
	if (k > 1) {
#line 681 "clasyf_rk.f"
	    i__1 = *n - k + 1;
#line 681 "clasyf_rk.f"
	    i__2 = k - 1;
#line 681 "clasyf_rk.f"
	    z__1.r = -1., z__1.i = -0.;
#line 681 "clasyf_rk.f"
	    cgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &
		    w[k + w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1, (
		    ftnlen)12);
#line 681 "clasyf_rk.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 688 "clasyf_rk.f"
	i__1 = k + k * w_dim1;
#line 688 "clasyf_rk.f"
	absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[k + k * 
		w_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 694 "clasyf_rk.f"
	if (k < *n) {
#line 695 "clasyf_rk.f"
	    i__1 = *n - k;
#line 695 "clasyf_rk.f"
	    imax = k + icamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 696 "clasyf_rk.f"
	    i__1 = imax + k * w_dim1;
#line 696 "clasyf_rk.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    k * w_dim1]), abs(d__2));
#line 697 "clasyf_rk.f"
	} else {
#line 698 "clasyf_rk.f"
	    colmax = 0.;
#line 699 "clasyf_rk.f"
	}

#line 701 "clasyf_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 705 "clasyf_rk.f"
	    if (*info == 0) {
#line 705 "clasyf_rk.f"
		*info = k;
#line 705 "clasyf_rk.f"
	    }
#line 707 "clasyf_rk.f"
	    kp = k;
#line 708 "clasyf_rk.f"
	    i__1 = *n - k + 1;
#line 708 "clasyf_rk.f"
	    ccopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
		    c__1);

/*           Set E( K ) to zero */

#line 712 "clasyf_rk.f"
	    if (k < *n) {
#line 712 "clasyf_rk.f"
		i__1 = k;
#line 712 "clasyf_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 712 "clasyf_rk.f"
	    }

#line 715 "clasyf_rk.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 724 "clasyf_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 728 "clasyf_rk.f"
		kp = k;

#line 730 "clasyf_rk.f"
	    } else {

#line 732 "clasyf_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 736 "clasyf_rk.f"
L72:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column K+1 of W and update it */

#line 743 "clasyf_rk.f"
		i__1 = imax - k;
#line 743 "clasyf_rk.f"
		ccopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 744 "clasyf_rk.f"
		i__1 = *n - imax + 1;
#line 744 "clasyf_rk.f"
		ccopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 
			1) * w_dim1], &c__1);
#line 746 "clasyf_rk.f"
		if (k > 1) {
#line 746 "clasyf_rk.f"
		    i__1 = *n - k + 1;
#line 746 "clasyf_rk.f"
		    i__2 = k - 1;
#line 746 "clasyf_rk.f"
		    z__1.r = -1., z__1.i = -0.;
#line 746 "clasyf_rk.f"
		    cgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1]
			    , lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 
			    1) * w_dim1], &c__1, (ftnlen)12);
#line 746 "clasyf_rk.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 755 "clasyf_rk.f"
		if (imax != k) {
#line 756 "clasyf_rk.f"
		    i__1 = imax - k;
#line 756 "clasyf_rk.f"
		    jmax = k - 1 + icamax_(&i__1, &w[k + (k + 1) * w_dim1], &
			    c__1);
#line 757 "clasyf_rk.f"
		    i__1 = jmax + (k + 1) * w_dim1;
#line 757 "clasyf_rk.f"
		    rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    w[jmax + (k + 1) * w_dim1]), abs(d__2));
#line 758 "clasyf_rk.f"
		} else {
#line 759 "clasyf_rk.f"
		    rowmax = 0.;
#line 760 "clasyf_rk.f"
		}

#line 762 "clasyf_rk.f"
		if (imax < *n) {
#line 763 "clasyf_rk.f"
		    i__1 = *n - imax;
#line 763 "clasyf_rk.f"
		    itemp = imax + icamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
#line 764 "clasyf_rk.f"
		    i__1 = itemp + (k + 1) * w_dim1;
#line 764 "clasyf_rk.f"
		    stemp = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    itemp + (k + 1) * w_dim1]), abs(d__2));
#line 765 "clasyf_rk.f"
		    if (stemp > rowmax) {
#line 766 "clasyf_rk.f"
			rowmax = stemp;
#line 767 "clasyf_rk.f"
			jmax = itemp;
#line 768 "clasyf_rk.f"
		    }
#line 769 "clasyf_rk.f"
		}

/*                 Equivalent to testing for */
/*                 CABS1( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 775 "clasyf_rk.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 775 "clasyf_rk.f"
		if (! ((d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax 
			+ (k + 1) * w_dim1]), abs(d__2)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 781 "clasyf_rk.f"
		    kp = imax;

/*                    copy column K+1 of W to column K of W */

#line 785 "clasyf_rk.f"
		    i__1 = *n - k + 1;
#line 785 "clasyf_rk.f"
		    ccopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 787 "clasyf_rk.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 792 "clasyf_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 798 "clasyf_rk.f"
		    kp = imax;
#line 799 "clasyf_rk.f"
		    kstep = 2;
#line 800 "clasyf_rk.f"
		    done = TRUE_;
#line 801 "clasyf_rk.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 805 "clasyf_rk.f"
		    p = imax;
#line 806 "clasyf_rk.f"
		    colmax = rowmax;
#line 807 "clasyf_rk.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 811 "clasyf_rk.f"
		    i__1 = *n - k + 1;
#line 811 "clasyf_rk.f"
		    ccopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 813 "clasyf_rk.f"
		}

/*                 End pivot search loop body */

#line 817 "clasyf_rk.f"
		if (! done) {
#line 817 "clasyf_rk.f"
		    goto L72;
#line 817 "clasyf_rk.f"
		}

#line 819 "clasyf_rk.f"
	    }

/*           ============================================================ */

#line 823 "clasyf_rk.f"
	    kk = k + kstep - 1;

#line 825 "clasyf_rk.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

#line 829 "clasyf_rk.f"
		i__1 = p - k;
#line 829 "clasyf_rk.f"
		ccopy_(&i__1, &a[k + k * a_dim1], &c__1, &a[p + k * a_dim1], 
			lda);
#line 830 "clasyf_rk.f"
		i__1 = *n - p + 1;
#line 830 "clasyf_rk.f"
		ccopy_(&i__1, &a[p + k * a_dim1], &c__1, &a[p + p * a_dim1], &
			c__1);

/*              Interchange rows K and P in first K columns of A */
/*              and first K+1 columns of W */

#line 835 "clasyf_rk.f"
		cswap_(&k, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 836 "clasyf_rk.f"
		cswap_(&kk, &w[k + w_dim1], ldw, &w[p + w_dim1], ldw);
#line 837 "clasyf_rk.f"
	    }

/*           Updated column KP is already stored in column KK of W */

#line 841 "clasyf_rk.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

#line 845 "clasyf_rk.f"
		i__1 = kp + k * a_dim1;
#line 845 "clasyf_rk.f"
		i__2 = kk + k * a_dim1;
#line 845 "clasyf_rk.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 846 "clasyf_rk.f"
		i__1 = kp - k - 1;
#line 846 "clasyf_rk.f"
		ccopy_(&i__1, &a[k + 1 + kk * a_dim1], &c__1, &a[kp + (k + 1) 
			* a_dim1], lda);
#line 847 "clasyf_rk.f"
		i__1 = *n - kp + 1;
#line 847 "clasyf_rk.f"
		ccopy_(&i__1, &a[kp + kk * a_dim1], &c__1, &a[kp + kp * 
			a_dim1], &c__1);

/*              Interchange rows KK and KP in first KK columns of A and W */

#line 851 "clasyf_rk.f"
		cswap_(&kk, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 852 "clasyf_rk.f"
		cswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 853 "clasyf_rk.f"
	    }

#line 855 "clasyf_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of W now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

/*              Store L(k) in column k of A */

#line 865 "clasyf_rk.f"
		i__1 = *n - k + 1;
#line 865 "clasyf_rk.f"
		ccopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 866 "clasyf_rk.f"
		if (k < *n) {
#line 867 "clasyf_rk.f"
		    i__1 = k + k * a_dim1;
#line 867 "clasyf_rk.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + 
			    k * a_dim1]), abs(d__2)) >= sfmin) {
#line 868 "clasyf_rk.f"
			z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 868 "clasyf_rk.f"
			r1.r = z__1.r, r1.i = z__1.i;
#line 869 "clasyf_rk.f"
			i__1 = *n - k;
#line 869 "clasyf_rk.f"
			cscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 870 "clasyf_rk.f"
		    } else /* if(complicated condition) */ {
#line 870 "clasyf_rk.f"
			i__1 = k + k * a_dim1;
#line 870 "clasyf_rk.f"
			if (a[i__1].r != 0. || a[i__1].i != 0.) {
#line 871 "clasyf_rk.f"
			    i__1 = *n;
#line 871 "clasyf_rk.f"
			    for (ii = k + 1; ii <= i__1; ++ii) {
#line 872 "clasyf_rk.f"
				i__2 = ii + k * a_dim1;
#line 872 "clasyf_rk.f"
				z_div(&z__1, &a[ii + k * a_dim1], &a[k + k * 
					a_dim1]);
#line 872 "clasyf_rk.f"
				a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 873 "clasyf_rk.f"
/* L74: */
#line 873 "clasyf_rk.f"
			    }
#line 874 "clasyf_rk.f"
			}
#line 874 "clasyf_rk.f"
		    }

/*                 Store the subdiagonal element of D in array E */

#line 878 "clasyf_rk.f"
		    i__1 = k;
#line 878 "clasyf_rk.f"
		    e[i__1].r = 0., e[i__1].i = 0.;

#line 880 "clasyf_rk.f"
		}

#line 882 "clasyf_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

#line 891 "clasyf_rk.f"
		if (k < *n - 1) {

/*                 Store L(k) and L(k+1) in columns k and k+1 of A */

#line 895 "clasyf_rk.f"
		    i__1 = k + 1 + k * w_dim1;
#line 895 "clasyf_rk.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 896 "clasyf_rk.f"
		    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
#line 896 "clasyf_rk.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 897 "clasyf_rk.f"
		    z_div(&z__1, &w[k + k * w_dim1], &d21);
#line 897 "clasyf_rk.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 898 "clasyf_rk.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 898 "clasyf_rk.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 898 "clasyf_rk.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 898 "clasyf_rk.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 899 "clasyf_rk.f"
		    i__1 = *n;
#line 899 "clasyf_rk.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 900 "clasyf_rk.f"
			i__2 = j + k * a_dim1;
#line 900 "clasyf_rk.f"
			i__3 = j + k * w_dim1;
#line 900 "clasyf_rk.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 900 "clasyf_rk.f"
			i__4 = j + (k + 1) * w_dim1;
#line 900 "clasyf_rk.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 900 "clasyf_rk.f"
			z_div(&z__2, &z__3, &d21);
#line 900 "clasyf_rk.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 900 "clasyf_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 902 "clasyf_rk.f"
			i__2 = j + (k + 1) * a_dim1;
#line 902 "clasyf_rk.f"
			i__3 = j + (k + 1) * w_dim1;
#line 902 "clasyf_rk.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 902 "clasyf_rk.f"
			i__4 = j + k * w_dim1;
#line 902 "clasyf_rk.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 902 "clasyf_rk.f"
			z_div(&z__2, &z__3, &d21);
#line 902 "clasyf_rk.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 902 "clasyf_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 904 "clasyf_rk.f"
/* L80: */
#line 904 "clasyf_rk.f"
		    }
#line 905 "clasyf_rk.f"
		}

/*              Copy diagonal elements of D(K) to A, */
/*              copy subdiagonal element of D(K) to E(K) and */
/*              ZERO out subdiagonal entry of A */

#line 911 "clasyf_rk.f"
		i__1 = k + k * a_dim1;
#line 911 "clasyf_rk.f"
		i__2 = k + k * w_dim1;
#line 911 "clasyf_rk.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 912 "clasyf_rk.f"
		i__1 = k + 1 + k * a_dim1;
#line 912 "clasyf_rk.f"
		a[i__1].r = 0., a[i__1].i = 0.;
#line 913 "clasyf_rk.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 913 "clasyf_rk.f"
		i__2 = k + 1 + (k + 1) * w_dim1;
#line 913 "clasyf_rk.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 914 "clasyf_rk.f"
		i__1 = k;
#line 914 "clasyf_rk.f"
		i__2 = k + 1 + k * w_dim1;
#line 914 "clasyf_rk.f"
		e[i__1].r = w[i__2].r, e[i__1].i = w[i__2].i;
#line 915 "clasyf_rk.f"
		i__1 = k + 1;
#line 915 "clasyf_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;

#line 917 "clasyf_rk.f"
	    }

/*           End column K is nonsingular */

#line 921 "clasyf_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 925 "clasyf_rk.f"
	if (kstep == 1) {
#line 926 "clasyf_rk.f"
	    ipiv[k] = kp;
#line 927 "clasyf_rk.f"
	} else {
#line 928 "clasyf_rk.f"
	    ipiv[k] = -p;
#line 929 "clasyf_rk.f"
	    ipiv[k + 1] = -kp;
#line 930 "clasyf_rk.f"
	}

/*        Increase K and return to the start of the main loop */

#line 934 "clasyf_rk.f"
	k += kstep;
#line 935 "clasyf_rk.f"
	goto L70;

#line 937 "clasyf_rk.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T */

/*        computing blocks of NB columns at a time */

#line 945 "clasyf_rk.f"
	i__1 = *n;
#line 945 "clasyf_rk.f"
	i__2 = *nb;
#line 945 "clasyf_rk.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 946 "clasyf_rk.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 946 "clasyf_rk.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 950 "clasyf_rk.f"
	    i__3 = j + jb - 1;
#line 950 "clasyf_rk.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 951 "clasyf_rk.f"
		i__4 = j + jb - jj;
#line 951 "clasyf_rk.f"
		i__5 = k - 1;
#line 951 "clasyf_rk.f"
		z__1.r = -1., z__1.i = -0.;
#line 951 "clasyf_rk.f"
		cgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1]
			, &c__1, (ftnlen)12);
#line 954 "clasyf_rk.f"
/* L100: */
#line 954 "clasyf_rk.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 958 "clasyf_rk.f"
	    if (j + jb <= *n) {
#line 958 "clasyf_rk.f"
		i__3 = *n - j - jb + 1;
#line 958 "clasyf_rk.f"
		i__4 = k - 1;
#line 958 "clasyf_rk.f"
		z__1.r = -1., z__1.i = -0.;
#line 958 "clasyf_rk.f"
		cgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, 
			&a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 958 "clasyf_rk.f"
	    }
#line 962 "clasyf_rk.f"
/* L110: */
#line 962 "clasyf_rk.f"
	}

/*        Set KB to the number of columns factorized */

#line 966 "clasyf_rk.f"
	*kb = k - 1;

#line 968 "clasyf_rk.f"
    }

#line 970 "clasyf_rk.f"
    return 0;

/*     End of CLASYF_RK */

} /* clasyf_rk__ */

