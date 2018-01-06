#line 1 "zlahef_rk.f"
/* zlahef_rk.f -- translated by f2c (version 20100827).
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

#line 1 "zlahef_rk.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLAHEF_RK computes a partial factorization of a complex Hermitian indefinite matrix using bound
ed Bunch-Kaufman (rook) diagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAHEF_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahef_
rk.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahef_
rk.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahef_
rk.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAHEF_RK( UPLO, N, NB, KB, A, LDA, E, IPIV, W, LDW, */
/*                             INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KB, LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), E( * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > ZLAHEF_RK computes a partial factorization of a complex Hermitian */
/* > matrix A using the bounded Bunch-Kaufman (rook) diagonal */
/* > pivoting method. The partial factorization has the form: */
/* > */
/* > A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or: */
/* >       ( 0  U22 ) (  0   D  ) ( U12**H U22**H ) */
/* > */
/* > A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L', */
/* >       ( L21  I ) (  0  A22 ) (  0       I    ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > */
/* > ZLAHEF_RK is an auxiliary routine called by ZHETRF_RK. It uses */
/* > blocked code (calling Level 3 BLAS) to update the submatrix */
/* > A11 (if UPLO = 'U') or A22 (if UPLO = 'L'). */
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
/* >          W is COMPLEX*16 array, dimension (LDW,NB) */
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

/* > \ingroup complex16HEcomputational */

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
/* Subroutine */ int zlahef_rk__(char *uplo, integer *n, integer *nb, integer 
	*kb, doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	doublecomplex *w, integer *ldw, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, p;
    static doublereal t, r1;
    static doublecomplex d11, d21, d22;
    static integer jb, ii, jj, kk, kp, kw, kkw;
    static logical done;
    static integer imax, jmax;
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dtemp, sfmin;
    static integer itemp;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer kstep;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal absakk;
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal colmax;
    extern /* Subroutine */ int zlacgv_(integer *, doublecomplex *, integer *)
	    ;
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

#line 319 "zlahef_rk.f"
    /* Parameter adjustments */
#line 319 "zlahef_rk.f"
    a_dim1 = *lda;
#line 319 "zlahef_rk.f"
    a_offset = 1 + a_dim1;
#line 319 "zlahef_rk.f"
    a -= a_offset;
#line 319 "zlahef_rk.f"
    --e;
#line 319 "zlahef_rk.f"
    --ipiv;
#line 319 "zlahef_rk.f"
    w_dim1 = *ldw;
#line 319 "zlahef_rk.f"
    w_offset = 1 + w_dim1;
#line 319 "zlahef_rk.f"
    w -= w_offset;
#line 319 "zlahef_rk.f"

#line 319 "zlahef_rk.f"
    /* Function Body */
#line 319 "zlahef_rk.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 323 "zlahef_rk.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 327 "zlahef_rk.f"
    sfmin = dlamch_("S", (ftnlen)1);

#line 329 "zlahef_rk.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 (note that conjg(W) is actually stored) */
/*        Initilize the first entry of array E, where superdiagonal */
/*        elements of D are stored */

#line 337 "zlahef_rk.f"
	e[1].r = 0., e[1].i = 0.;

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

#line 341 "zlahef_rk.f"
	k = *n;
#line 342 "zlahef_rk.f"
L10:

/*        KW is the column of W which corresponds to column K of A */

#line 346 "zlahef_rk.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 350 "zlahef_rk.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 350 "zlahef_rk.f"
	    goto L30;
#line 350 "zlahef_rk.f"
	}

#line 353 "zlahef_rk.f"
	kstep = 1;
#line 354 "zlahef_rk.f"
	p = k;

/*        Copy column K of A to column KW of W and update it */

#line 358 "zlahef_rk.f"
	if (k > 1) {
#line 358 "zlahef_rk.f"
	    i__1 = k - 1;
#line 358 "zlahef_rk.f"
	    zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &
		    c__1);
#line 358 "zlahef_rk.f"
	}
#line 360 "zlahef_rk.f"
	i__1 = k + kw * w_dim1;
#line 360 "zlahef_rk.f"
	i__2 = k + k * a_dim1;
#line 360 "zlahef_rk.f"
	d__1 = a[i__2].r;
#line 360 "zlahef_rk.f"
	w[i__1].r = d__1, w[i__1].i = 0.;
#line 361 "zlahef_rk.f"
	if (k < *n) {
#line 362 "zlahef_rk.f"
	    i__1 = *n - k;
#line 362 "zlahef_rk.f"
	    z__1.r = -1., z__1.i = -0.;
#line 362 "zlahef_rk.f"
	    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 364 "zlahef_rk.f"
	    i__1 = k + kw * w_dim1;
#line 364 "zlahef_rk.f"
	    i__2 = k + kw * w_dim1;
#line 364 "zlahef_rk.f"
	    d__1 = w[i__2].r;
#line 364 "zlahef_rk.f"
	    w[i__1].r = d__1, w[i__1].i = 0.;
#line 365 "zlahef_rk.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 370 "zlahef_rk.f"
	i__1 = k + kw * w_dim1;
#line 370 "zlahef_rk.f"
	absakk = (d__1 = w[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 376 "zlahef_rk.f"
	if (k > 1) {
#line 377 "zlahef_rk.f"
	    i__1 = k - 1;
#line 377 "zlahef_rk.f"
	    imax = izamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 378 "zlahef_rk.f"
	    i__1 = imax + kw * w_dim1;
#line 378 "zlahef_rk.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    kw * w_dim1]), abs(d__2));
#line 379 "zlahef_rk.f"
	} else {
#line 380 "zlahef_rk.f"
	    colmax = 0.;
#line 381 "zlahef_rk.f"
	}

#line 383 "zlahef_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 387 "zlahef_rk.f"
	    if (*info == 0) {
#line 387 "zlahef_rk.f"
		*info = k;
#line 387 "zlahef_rk.f"
	    }
#line 389 "zlahef_rk.f"
	    kp = k;
#line 390 "zlahef_rk.f"
	    i__1 = k + k * a_dim1;
#line 390 "zlahef_rk.f"
	    i__2 = k + kw * w_dim1;
#line 390 "zlahef_rk.f"
	    d__1 = w[i__2].r;
#line 390 "zlahef_rk.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 391 "zlahef_rk.f"
	    if (k > 1) {
#line 391 "zlahef_rk.f"
		i__1 = k - 1;
#line 391 "zlahef_rk.f"
		zcopy_(&i__1, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], 
			&c__1);
#line 391 "zlahef_rk.f"
	    }

/*           Set E( K ) to zero */

#line 396 "zlahef_rk.f"
	    if (k > 1) {
#line 396 "zlahef_rk.f"
		i__1 = k;
#line 396 "zlahef_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 396 "zlahef_rk.f"
	    }

#line 399 "zlahef_rk.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */
#line 408 "zlahef_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 412 "zlahef_rk.f"
		kp = k;

#line 414 "zlahef_rk.f"
	    } else {

/*              Lop until pivot found */

#line 418 "zlahef_rk.f"
		done = FALSE_;

#line 420 "zlahef_rk.f"
L12:

/*                 BEGIN pivot search loop body */


/*                 Copy column IMAX to column KW-1 of W and update it */

#line 427 "zlahef_rk.f"
		if (imax > 1) {
#line 427 "zlahef_rk.f"
		    i__1 = imax - 1;
#line 427 "zlahef_rk.f"
		    zcopy_(&i__1, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			    w_dim1 + 1], &c__1);
#line 427 "zlahef_rk.f"
		}
#line 430 "zlahef_rk.f"
		i__1 = imax + (kw - 1) * w_dim1;
#line 430 "zlahef_rk.f"
		i__2 = imax + imax * a_dim1;
#line 430 "zlahef_rk.f"
		d__1 = a[i__2].r;
#line 430 "zlahef_rk.f"
		w[i__1].r = d__1, w[i__1].i = 0.;

#line 432 "zlahef_rk.f"
		i__1 = k - imax;
#line 432 "zlahef_rk.f"
		zcopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);
#line 434 "zlahef_rk.f"
		i__1 = k - imax;
#line 434 "zlahef_rk.f"
		zlacgv_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);

#line 436 "zlahef_rk.f"
		if (k < *n) {
#line 437 "zlahef_rk.f"
		    i__1 = *n - k;
#line 437 "zlahef_rk.f"
		    z__1.r = -1., z__1.i = -0.;
#line 437 "zlahef_rk.f"
		    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 440 "zlahef_rk.f"
		    i__1 = imax + (kw - 1) * w_dim1;
#line 440 "zlahef_rk.f"
		    i__2 = imax + (kw - 1) * w_dim1;
#line 440 "zlahef_rk.f"
		    d__1 = w[i__2].r;
#line 440 "zlahef_rk.f"
		    w[i__1].r = d__1, w[i__1].i = 0.;
#line 441 "zlahef_rk.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 447 "zlahef_rk.f"
		if (imax != k) {
#line 448 "zlahef_rk.f"
		    i__1 = k - imax;
#line 448 "zlahef_rk.f"
		    jmax = imax + izamax_(&i__1, &w[imax + 1 + (kw - 1) * 
			    w_dim1], &c__1);
#line 450 "zlahef_rk.f"
		    i__1 = jmax + (kw - 1) * w_dim1;
#line 450 "zlahef_rk.f"
		    rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    w[jmax + (kw - 1) * w_dim1]), abs(d__2));
#line 451 "zlahef_rk.f"
		} else {
#line 452 "zlahef_rk.f"
		    rowmax = 0.;
#line 453 "zlahef_rk.f"
		}

#line 455 "zlahef_rk.f"
		if (imax > 1) {
#line 456 "zlahef_rk.f"
		    i__1 = imax - 1;
#line 456 "zlahef_rk.f"
		    itemp = izamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
#line 457 "zlahef_rk.f"
		    i__1 = itemp + (kw - 1) * w_dim1;
#line 457 "zlahef_rk.f"
		    dtemp = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    itemp + (kw - 1) * w_dim1]), abs(d__2));
#line 458 "zlahef_rk.f"
		    if (dtemp > rowmax) {
#line 459 "zlahef_rk.f"
			rowmax = dtemp;
#line 460 "zlahef_rk.f"
			jmax = itemp;
#line 461 "zlahef_rk.f"
		    }
#line 462 "zlahef_rk.f"
		}

/*                 Case(2) */
/*                 Equivalent to testing for */
/*                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 469 "zlahef_rk.f"
		i__1 = imax + (kw - 1) * w_dim1;
#line 469 "zlahef_rk.f"
		if (! ((d__1 = w[i__1].r, abs(d__1)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 475 "zlahef_rk.f"
		    kp = imax;

/*                    copy column KW-1 of W to column KW of W */

#line 479 "zlahef_rk.f"
		    zcopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 481 "zlahef_rk.f"
		    done = TRUE_;

/*                 Case(3) */
/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 487 "zlahef_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 493 "zlahef_rk.f"
		    kp = imax;
#line 494 "zlahef_rk.f"
		    kstep = 2;
#line 495 "zlahef_rk.f"
		    done = TRUE_;

/*                 Case(4) */
#line 498 "zlahef_rk.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 502 "zlahef_rk.f"
		    p = imax;
#line 503 "zlahef_rk.f"
		    colmax = rowmax;
#line 504 "zlahef_rk.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 508 "zlahef_rk.f"
		    zcopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 510 "zlahef_rk.f"
		}


/*                 END pivot search loop body */

#line 515 "zlahef_rk.f"
		if (! done) {
#line 515 "zlahef_rk.f"
		    goto L12;
#line 515 "zlahef_rk.f"
		}

#line 517 "zlahef_rk.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 525 "zlahef_rk.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 529 "zlahef_rk.f"
	    kkw = *nb + kk - *n;

/*           Interchange rows and columns P and K. */
/*           Updated column P is already stored in column KW of W. */

#line 534 "zlahef_rk.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P of submatrix A */
/*              at step K. No need to copy element into columns */
/*              K and K-1 of A for 2-by-2 pivot, since these columns */
/*              will be later overwritten. */

#line 541 "zlahef_rk.f"
		i__1 = p + p * a_dim1;
#line 541 "zlahef_rk.f"
		i__2 = k + k * a_dim1;
#line 541 "zlahef_rk.f"
		d__1 = a[i__2].r;
#line 541 "zlahef_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 542 "zlahef_rk.f"
		i__1 = k - 1 - p;
#line 542 "zlahef_rk.f"
		zcopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * 
			a_dim1], lda);
#line 544 "zlahef_rk.f"
		i__1 = k - 1 - p;
#line 544 "zlahef_rk.f"
		zlacgv_(&i__1, &a[p + (p + 1) * a_dim1], lda);
#line 545 "zlahef_rk.f"
		if (p > 1) {
#line 545 "zlahef_rk.f"
		    i__1 = p - 1;
#line 545 "zlahef_rk.f"
		    zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 
			    1], &c__1);
#line 545 "zlahef_rk.f"
		}

/*              Interchange rows K and P in the last K+1 to N columns of A */
/*              (columns K and K-1 of A for 2-by-2 pivot will be */
/*              later overwritten). Interchange rows K and P */
/*              in last KKW to NB columns of W. */

#line 553 "zlahef_rk.f"
		if (k < *n) {
#line 553 "zlahef_rk.f"
		    i__1 = *n - k;
#line 553 "zlahef_rk.f"
		    zswap_(&i__1, &a[k + (k + 1) * a_dim1], lda, &a[p + (k + 
			    1) * a_dim1], lda);
#line 553 "zlahef_rk.f"
		}
#line 556 "zlahef_rk.f"
		i__1 = *n - kk + 1;
#line 556 "zlahef_rk.f"
		zswap_(&i__1, &w[k + kkw * w_dim1], ldw, &w[p + kkw * w_dim1],
			 ldw);
#line 558 "zlahef_rk.f"
	    }

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KKW of W. */

#line 563 "zlahef_rk.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K-1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 570 "zlahef_rk.f"
		i__1 = kp + kp * a_dim1;
#line 570 "zlahef_rk.f"
		i__2 = kk + kk * a_dim1;
#line 570 "zlahef_rk.f"
		d__1 = a[i__2].r;
#line 570 "zlahef_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 571 "zlahef_rk.f"
		i__1 = kk - 1 - kp;
#line 571 "zlahef_rk.f"
		zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 573 "zlahef_rk.f"
		i__1 = kk - 1 - kp;
#line 573 "zlahef_rk.f"
		zlacgv_(&i__1, &a[kp + (kp + 1) * a_dim1], lda);
#line 574 "zlahef_rk.f"
		if (kp > 1) {
#line 574 "zlahef_rk.f"
		    i__1 = kp - 1;
#line 574 "zlahef_rk.f"
		    zcopy_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 574 "zlahef_rk.f"
		}

/*              Interchange rows KK and KP in last K+1 to N columns of A */
/*              (columns K (or K and K-1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in last KKW to NB columns of W. */

#line 582 "zlahef_rk.f"
		if (k < *n) {
#line 582 "zlahef_rk.f"
		    i__1 = *n - k;
#line 582 "zlahef_rk.f"
		    zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 582 "zlahef_rk.f"
		}
#line 585 "zlahef_rk.f"
		i__1 = *n - kk + 1;
#line 585 "zlahef_rk.f"
		zswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 587 "zlahef_rk.f"
	    }

#line 589 "zlahef_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column kw of W now holds */

/*              W(kw) = U(k)*D(k), */

/*              where U(k) is the k-th column of U */

/*              (1) Store subdiag. elements of column U(k) */
/*              and 1-by-1 block D(k) in column k of A. */
/*              (NOTE: Diagonal element U(k,k) is a UNIT element */
/*              and not stored) */
/*                 A(k,k) := D(k,k) = W(k,kw) */
/*                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k) */

/*              (NOTE: No need to use for Hermitian matrix */
/*              A( K, K ) = REAL( W( K, K) ) to separately copy diagonal */
/*              element D(k,k) from W (potentially saves only one load)) */
#line 607 "zlahef_rk.f"
		zcopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 608 "zlahef_rk.f"
		if (k > 1) {

/*                 (NOTE: No need to check if A(k,k) is NOT ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  case A(k,k) = 0 falls into 2x2 pivot case(3)) */

/*                 Handle division by a small number */

#line 616 "zlahef_rk.f"
		    i__1 = k + k * a_dim1;
#line 616 "zlahef_rk.f"
		    t = a[i__1].r;
#line 617 "zlahef_rk.f"
		    if (abs(t) >= sfmin) {
#line 618 "zlahef_rk.f"
			r1 = 1. / t;
#line 619 "zlahef_rk.f"
			i__1 = k - 1;
#line 619 "zlahef_rk.f"
			zdscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 620 "zlahef_rk.f"
		    } else {
#line 621 "zlahef_rk.f"
			i__1 = k - 1;
#line 621 "zlahef_rk.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 622 "zlahef_rk.f"
			    i__2 = ii + k * a_dim1;
#line 622 "zlahef_rk.f"
			    i__3 = ii + k * a_dim1;
#line 622 "zlahef_rk.f"
			    z__1.r = a[i__3].r / t, z__1.i = a[i__3].i / t;
#line 622 "zlahef_rk.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 623 "zlahef_rk.f"
/* L14: */
#line 623 "zlahef_rk.f"
			}
#line 624 "zlahef_rk.f"
		    }

/*                 (2) Conjugate column W(kw) */

#line 628 "zlahef_rk.f"
		    i__1 = k - 1;
#line 628 "zlahef_rk.f"
		    zlacgv_(&i__1, &w[kw * w_dim1 + 1], &c__1);

/*                 Store the superdiagonal element of D in array E */

#line 632 "zlahef_rk.f"
		    i__1 = k;
#line 632 "zlahef_rk.f"
		    e[i__1].r = 0., e[i__1].i = 0.;

#line 634 "zlahef_rk.f"
		}

#line 636 "zlahef_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns kw and kw-1 of W now hold */

/*              ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              (1) Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2 */
/*              block D(k-1:k,k-1:k) in columns k-1 and k of A. */
/*              (NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT */
/*              block and not stored) */
/*                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw) */
/*                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) = */
/*                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) ) */

#line 653 "zlahef_rk.f"
		if (k > 2) {

/*                 Factor out the columns of the inverse of 2-by-2 pivot */
/*                 block D, so that each column contains 1, to reduce the */
/*                 number of FLOPS when we multiply panel */
/*                 ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1). */

/*                 D**(-1) = ( d11 cj(d21) )**(-1) = */
/*                           ( d21    d22 ) */

/*                 = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) = */
/*                                          ( (-d21) (     d11 ) ) */

/*                 = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) * */

/*                   * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) = */
/*                     (     (      -1 )           ( d11/conj(d21) ) ) */

/*                 = 1/(|d21|**2) * 1/(D22*D11-1) * */

/*                   * ( d21*( D11 ) conj(d21)*(  -1 ) ) = */
/*                     (     (  -1 )           ( D22 ) ) */

/*                 = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) = */
/*                                      (     (  -1 )           ( D22 ) ) */

/*                 = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) = */
/*                   (               (  -1 )         ( D22 ) ) */

/*                 Handle division by a small number. (NOTE: order of */
/*                 operations is important) */

/*                 = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) ) */
/*                   (   ((  -1 )          )   (( D22 )     ) ), */

/*                 where D11 = d22/d21, */
/*                       D22 = d11/conj(d21), */
/*                       D21 = d21, */
/*                       T = 1/(D22*D11-1). */

/*                 (NOTE: No need to check for division by ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  (a) d21 != 0 in 2x2 pivot case(4), */
/*                      since |d21| should be larger than |d11| and |d22|; */
/*                  (b) (D22*D11 - 1) != 0, since from (a), */
/*                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.) */

#line 700 "zlahef_rk.f"
		    i__1 = k - 1 + kw * w_dim1;
#line 700 "zlahef_rk.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 701 "zlahef_rk.f"
		    d_cnjg(&z__2, &d21);
#line 701 "zlahef_rk.f"
		    z_div(&z__1, &w[k + kw * w_dim1], &z__2);
#line 701 "zlahef_rk.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 702 "zlahef_rk.f"
		    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d21);
#line 702 "zlahef_rk.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 703 "zlahef_rk.f"
		    z__1.r = d11.r * d22.r - d11.i * d22.i, z__1.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 703 "zlahef_rk.f"
		    t = 1. / (z__1.r - 1.);

/*                 Update elements in columns A(k-1) and A(k) as */
/*                 dot products of rows of ( W(kw-1) W(kw) ) and columns */
/*                 of D**(-1) */

#line 709 "zlahef_rk.f"
		    i__1 = k - 2;
#line 709 "zlahef_rk.f"
		    for (j = 1; j <= i__1; ++j) {
#line 710 "zlahef_rk.f"
			i__2 = j + (k - 1) * a_dim1;
#line 710 "zlahef_rk.f"
			i__3 = j + (kw - 1) * w_dim1;
#line 710 "zlahef_rk.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 710 "zlahef_rk.f"
			i__4 = j + kw * w_dim1;
#line 710 "zlahef_rk.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 710 "zlahef_rk.f"
			z_div(&z__2, &z__3, &d21);
#line 710 "zlahef_rk.f"
			z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 710 "zlahef_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 712 "zlahef_rk.f"
			i__2 = j + k * a_dim1;
#line 712 "zlahef_rk.f"
			i__3 = j + kw * w_dim1;
#line 712 "zlahef_rk.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 712 "zlahef_rk.f"
			i__4 = j + (kw - 1) * w_dim1;
#line 712 "zlahef_rk.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 712 "zlahef_rk.f"
			d_cnjg(&z__5, &d21);
#line 712 "zlahef_rk.f"
			z_div(&z__2, &z__3, &z__5);
#line 712 "zlahef_rk.f"
			z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 712 "zlahef_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 714 "zlahef_rk.f"
/* L20: */
#line 714 "zlahef_rk.f"
		    }
#line 715 "zlahef_rk.f"
		}

/*              Copy diagonal elements of D(K) to A, */
/*              copy superdiagonal element of D(K) to E(K) and */
/*              ZERO out superdiagonal entry of A */

#line 721 "zlahef_rk.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 721 "zlahef_rk.f"
		i__2 = k - 1 + (kw - 1) * w_dim1;
#line 721 "zlahef_rk.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 722 "zlahef_rk.f"
		i__1 = k - 1 + k * a_dim1;
#line 722 "zlahef_rk.f"
		a[i__1].r = 0., a[i__1].i = 0.;
#line 723 "zlahef_rk.f"
		i__1 = k + k * a_dim1;
#line 723 "zlahef_rk.f"
		i__2 = k + kw * w_dim1;
#line 723 "zlahef_rk.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 724 "zlahef_rk.f"
		i__1 = k;
#line 724 "zlahef_rk.f"
		i__2 = k - 1 + kw * w_dim1;
#line 724 "zlahef_rk.f"
		e[i__1].r = w[i__2].r, e[i__1].i = w[i__2].i;
#line 725 "zlahef_rk.f"
		i__1 = k - 1;
#line 725 "zlahef_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;

/*              (2) Conjugate columns W(kw) and W(kw-1) */

#line 729 "zlahef_rk.f"
		i__1 = k - 1;
#line 729 "zlahef_rk.f"
		zlacgv_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 730 "zlahef_rk.f"
		i__1 = k - 2;
#line 730 "zlahef_rk.f"
		zlacgv_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);

#line 732 "zlahef_rk.f"
	    }

/*           End column K is nonsingular */

#line 736 "zlahef_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 740 "zlahef_rk.f"
	if (kstep == 1) {
#line 741 "zlahef_rk.f"
	    ipiv[k] = kp;
#line 742 "zlahef_rk.f"
	} else {
#line 743 "zlahef_rk.f"
	    ipiv[k] = -p;
#line 744 "zlahef_rk.f"
	    ipiv[k - 1] = -kp;
#line 745 "zlahef_rk.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 749 "zlahef_rk.f"
	k -= kstep;
#line 750 "zlahef_rk.f"
	goto L10;

#line 752 "zlahef_rk.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**H = A11 - U12*W**H */

/*        computing blocks of NB columns at a time (note that conjg(W) is */
/*        actually stored) */

#line 761 "zlahef_rk.f"
	i__1 = -(*nb);
#line 761 "zlahef_rk.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 762 "zlahef_rk.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 762 "zlahef_rk.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 766 "zlahef_rk.f"
	    i__2 = j + jb - 1;
#line 766 "zlahef_rk.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 767 "zlahef_rk.f"
		i__3 = jj + jj * a_dim1;
#line 767 "zlahef_rk.f"
		i__4 = jj + jj * a_dim1;
#line 767 "zlahef_rk.f"
		d__1 = a[i__4].r;
#line 767 "zlahef_rk.f"
		a[i__3].r = d__1, a[i__3].i = 0.;
#line 768 "zlahef_rk.f"
		i__3 = jj - j + 1;
#line 768 "zlahef_rk.f"
		i__4 = *n - k;
#line 768 "zlahef_rk.f"
		z__1.r = -1., z__1.i = -0.;
#line 768 "zlahef_rk.f"
		zgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, 
			&a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 771 "zlahef_rk.f"
		i__3 = jj + jj * a_dim1;
#line 771 "zlahef_rk.f"
		i__4 = jj + jj * a_dim1;
#line 771 "zlahef_rk.f"
		d__1 = a[i__4].r;
#line 771 "zlahef_rk.f"
		a[i__3].r = d__1, a[i__3].i = 0.;
#line 772 "zlahef_rk.f"
/* L40: */
#line 772 "zlahef_rk.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 776 "zlahef_rk.f"
	    if (j >= 2) {
#line 776 "zlahef_rk.f"
		i__2 = j - 1;
#line 776 "zlahef_rk.f"
		i__3 = *n - k;
#line 776 "zlahef_rk.f"
		z__1.r = -1., z__1.i = -0.;
#line 776 "zlahef_rk.f"
		zgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, 
			&a[(k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * 
			w_dim1], ldw, &c_b1, &a[j * a_dim1 + 1], lda, (ftnlen)
			12, (ftnlen)9);
#line 776 "zlahef_rk.f"
	    }
#line 780 "zlahef_rk.f"
/* L50: */
#line 780 "zlahef_rk.f"
	}

/*        Set KB to the number of columns factorized */

#line 784 "zlahef_rk.f"
	*kb = *n - k;

#line 786 "zlahef_rk.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 (note that conjg(W) is actually stored) */

/*        Initilize the unused last entry of the subdiagonal array E. */

#line 794 "zlahef_rk.f"
	i__1 = *n;
#line 794 "zlahef_rk.f"
	e[i__1].r = 0., e[i__1].i = 0.;

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 798 "zlahef_rk.f"
	k = 1;
#line 799 "zlahef_rk.f"
L70:

/*        Exit from loop */

#line 803 "zlahef_rk.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 803 "zlahef_rk.f"
	    goto L90;
#line 803 "zlahef_rk.f"
	}

#line 806 "zlahef_rk.f"
	kstep = 1;
#line 807 "zlahef_rk.f"
	p = k;

/*        Copy column K of A to column K of W and update column K of W */

#line 811 "zlahef_rk.f"
	i__1 = k + k * w_dim1;
#line 811 "zlahef_rk.f"
	i__2 = k + k * a_dim1;
#line 811 "zlahef_rk.f"
	d__1 = a[i__2].r;
#line 811 "zlahef_rk.f"
	w[i__1].r = d__1, w[i__1].i = 0.;
#line 812 "zlahef_rk.f"
	if (k < *n) {
#line 812 "zlahef_rk.f"
	    i__1 = *n - k;
#line 812 "zlahef_rk.f"
	    zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &w[k + 1 + k * 
		    w_dim1], &c__1);
#line 812 "zlahef_rk.f"
	}
#line 814 "zlahef_rk.f"
	if (k > 1) {
#line 815 "zlahef_rk.f"
	    i__1 = *n - k + 1;
#line 815 "zlahef_rk.f"
	    i__2 = k - 1;
#line 815 "zlahef_rk.f"
	    z__1.r = -1., z__1.i = -0.;
#line 815 "zlahef_rk.f"
	    zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &
		    w[k + w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1, (
		    ftnlen)12);
#line 817 "zlahef_rk.f"
	    i__1 = k + k * w_dim1;
#line 817 "zlahef_rk.f"
	    i__2 = k + k * w_dim1;
#line 817 "zlahef_rk.f"
	    d__1 = w[i__2].r;
#line 817 "zlahef_rk.f"
	    w[i__1].r = d__1, w[i__1].i = 0.;
#line 818 "zlahef_rk.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 823 "zlahef_rk.f"
	i__1 = k + k * w_dim1;
#line 823 "zlahef_rk.f"
	absakk = (d__1 = w[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 829 "zlahef_rk.f"
	if (k < *n) {
#line 830 "zlahef_rk.f"
	    i__1 = *n - k;
#line 830 "zlahef_rk.f"
	    imax = k + izamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 831 "zlahef_rk.f"
	    i__1 = imax + k * w_dim1;
#line 831 "zlahef_rk.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    k * w_dim1]), abs(d__2));
#line 832 "zlahef_rk.f"
	} else {
#line 833 "zlahef_rk.f"
	    colmax = 0.;
#line 834 "zlahef_rk.f"
	}

#line 836 "zlahef_rk.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 840 "zlahef_rk.f"
	    if (*info == 0) {
#line 840 "zlahef_rk.f"
		*info = k;
#line 840 "zlahef_rk.f"
	    }
#line 842 "zlahef_rk.f"
	    kp = k;
#line 843 "zlahef_rk.f"
	    i__1 = k + k * a_dim1;
#line 843 "zlahef_rk.f"
	    i__2 = k + k * w_dim1;
#line 843 "zlahef_rk.f"
	    d__1 = w[i__2].r;
#line 843 "zlahef_rk.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 844 "zlahef_rk.f"
	    if (k < *n) {
#line 844 "zlahef_rk.f"
		i__1 = *n - k;
#line 844 "zlahef_rk.f"
		zcopy_(&i__1, &w[k + 1 + k * w_dim1], &c__1, &a[k + 1 + k * 
			a_dim1], &c__1);
#line 844 "zlahef_rk.f"
	    }

/*           Set E( K ) to zero */

#line 849 "zlahef_rk.f"
	    if (k < *n) {
#line 849 "zlahef_rk.f"
		i__1 = k;
#line 849 "zlahef_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;
#line 849 "zlahef_rk.f"
	    }

#line 852 "zlahef_rk.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 862 "zlahef_rk.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 866 "zlahef_rk.f"
		kp = k;

#line 868 "zlahef_rk.f"
	    } else {

#line 870 "zlahef_rk.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 874 "zlahef_rk.f"
L72:

/*                 BEGIN pivot search loop body */


/*                 Copy column IMAX to column k+1 of W and update it */

#line 881 "zlahef_rk.f"
		i__1 = imax - k;
#line 881 "zlahef_rk.f"
		zcopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 882 "zlahef_rk.f"
		i__1 = imax - k;
#line 882 "zlahef_rk.f"
		zlacgv_(&i__1, &w[k + (k + 1) * w_dim1], &c__1);
#line 883 "zlahef_rk.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 883 "zlahef_rk.f"
		i__2 = imax + imax * a_dim1;
#line 883 "zlahef_rk.f"
		d__1 = a[i__2].r;
#line 883 "zlahef_rk.f"
		w[i__1].r = d__1, w[i__1].i = 0.;

#line 885 "zlahef_rk.f"
		if (imax < *n) {
#line 885 "zlahef_rk.f"
		    i__1 = *n - imax;
#line 885 "zlahef_rk.f"
		    zcopy_(&i__1, &a[imax + 1 + imax * a_dim1], &c__1, &w[
			    imax + 1 + (k + 1) * w_dim1], &c__1);
#line 885 "zlahef_rk.f"
		}

#line 889 "zlahef_rk.f"
		if (k > 1) {
#line 890 "zlahef_rk.f"
		    i__1 = *n - k + 1;
#line 890 "zlahef_rk.f"
		    i__2 = k - 1;
#line 890 "zlahef_rk.f"
		    z__1.r = -1., z__1.i = -0.;
#line 890 "zlahef_rk.f"
		    zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1]
			    , lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 
			    1) * w_dim1], &c__1, (ftnlen)12);
#line 893 "zlahef_rk.f"
		    i__1 = imax + (k + 1) * w_dim1;
#line 893 "zlahef_rk.f"
		    i__2 = imax + (k + 1) * w_dim1;
#line 893 "zlahef_rk.f"
		    d__1 = w[i__2].r;
#line 893 "zlahef_rk.f"
		    w[i__1].r = d__1, w[i__1].i = 0.;
#line 894 "zlahef_rk.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 900 "zlahef_rk.f"
		if (imax != k) {
#line 901 "zlahef_rk.f"
		    i__1 = imax - k;
#line 901 "zlahef_rk.f"
		    jmax = k - 1 + izamax_(&i__1, &w[k + (k + 1) * w_dim1], &
			    c__1);
#line 902 "zlahef_rk.f"
		    i__1 = jmax + (k + 1) * w_dim1;
#line 902 "zlahef_rk.f"
		    rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    w[jmax + (k + 1) * w_dim1]), abs(d__2));
#line 903 "zlahef_rk.f"
		} else {
#line 904 "zlahef_rk.f"
		    rowmax = 0.;
#line 905 "zlahef_rk.f"
		}

#line 907 "zlahef_rk.f"
		if (imax < *n) {
#line 908 "zlahef_rk.f"
		    i__1 = *n - imax;
#line 908 "zlahef_rk.f"
		    itemp = imax + izamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
#line 909 "zlahef_rk.f"
		    i__1 = itemp + (k + 1) * w_dim1;
#line 909 "zlahef_rk.f"
		    dtemp = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    itemp + (k + 1) * w_dim1]), abs(d__2));
#line 910 "zlahef_rk.f"
		    if (dtemp > rowmax) {
#line 911 "zlahef_rk.f"
			rowmax = dtemp;
#line 912 "zlahef_rk.f"
			jmax = itemp;
#line 913 "zlahef_rk.f"
		    }
#line 914 "zlahef_rk.f"
		}

/*                 Case(2) */
/*                 Equivalent to testing for */
/*                 ABS( REAL( W( IMAX,K+1 ) ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 921 "zlahef_rk.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 921 "zlahef_rk.f"
		if (! ((d__1 = w[i__1].r, abs(d__1)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 927 "zlahef_rk.f"
		    kp = imax;

/*                    copy column K+1 of W to column K of W */

#line 931 "zlahef_rk.f"
		    i__1 = *n - k + 1;
#line 931 "zlahef_rk.f"
		    zcopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 933 "zlahef_rk.f"
		    done = TRUE_;

/*                 Case(3) */
/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 939 "zlahef_rk.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 945 "zlahef_rk.f"
		    kp = imax;
#line 946 "zlahef_rk.f"
		    kstep = 2;
#line 947 "zlahef_rk.f"
		    done = TRUE_;

/*                 Case(4) */
#line 950 "zlahef_rk.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 954 "zlahef_rk.f"
		    p = imax;
#line 955 "zlahef_rk.f"
		    colmax = rowmax;
#line 956 "zlahef_rk.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 960 "zlahef_rk.f"
		    i__1 = *n - k + 1;
#line 960 "zlahef_rk.f"
		    zcopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 962 "zlahef_rk.f"
		}


/*                 End pivot search loop body */

#line 967 "zlahef_rk.f"
		if (! done) {
#line 967 "zlahef_rk.f"
		    goto L72;
#line 967 "zlahef_rk.f"
		}

#line 969 "zlahef_rk.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 977 "zlahef_rk.f"
	    kk = k + kstep - 1;

/*           Interchange rows and columns P and K (only for 2-by-2 pivot). */
/*           Updated column P is already stored in column K of W. */

#line 982 "zlahef_rk.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column KK-1 to column P of submatrix A */
/*              at step K. No need to copy element into columns */
/*              K and K+1 of A for 2-by-2 pivot, since these columns */
/*              will be later overwritten. */

#line 989 "zlahef_rk.f"
		i__1 = p + p * a_dim1;
#line 989 "zlahef_rk.f"
		i__2 = k + k * a_dim1;
#line 989 "zlahef_rk.f"
		d__1 = a[i__2].r;
#line 989 "zlahef_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 990 "zlahef_rk.f"
		i__1 = p - k - 1;
#line 990 "zlahef_rk.f"
		zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[p + (k + 1) * 
			a_dim1], lda);
#line 991 "zlahef_rk.f"
		i__1 = p - k - 1;
#line 991 "zlahef_rk.f"
		zlacgv_(&i__1, &a[p + (k + 1) * a_dim1], lda);
#line 992 "zlahef_rk.f"
		if (p < *n) {
#line 992 "zlahef_rk.f"
		    i__1 = *n - p;
#line 992 "zlahef_rk.f"
		    zcopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p 
			    * a_dim1], &c__1);
#line 992 "zlahef_rk.f"
		}

/*              Interchange rows K and P in first K-1 columns of A */
/*              (columns K and K+1 of A for 2-by-2 pivot will be */
/*              later overwritten). Interchange rows K and P */
/*              in first KK columns of W. */

#line 1000 "zlahef_rk.f"
		if (k > 1) {
#line 1000 "zlahef_rk.f"
		    i__1 = k - 1;
#line 1000 "zlahef_rk.f"
		    zswap_(&i__1, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 1000 "zlahef_rk.f"
		}
#line 1002 "zlahef_rk.f"
		zswap_(&kk, &w[k + w_dim1], ldw, &w[p + w_dim1], ldw);
#line 1003 "zlahef_rk.f"
	    }

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KK of W. */

#line 1008 "zlahef_rk.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K+1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 1015 "zlahef_rk.f"
		i__1 = kp + kp * a_dim1;
#line 1015 "zlahef_rk.f"
		i__2 = kk + kk * a_dim1;
#line 1015 "zlahef_rk.f"
		d__1 = a[i__2].r;
#line 1015 "zlahef_rk.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 1016 "zlahef_rk.f"
		i__1 = kp - kk - 1;
#line 1016 "zlahef_rk.f"
		zcopy_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 1018 "zlahef_rk.f"
		i__1 = kp - kk - 1;
#line 1018 "zlahef_rk.f"
		zlacgv_(&i__1, &a[kp + (kk + 1) * a_dim1], lda);
#line 1019 "zlahef_rk.f"
		if (kp < *n) {
#line 1019 "zlahef_rk.f"
		    i__1 = *n - kp;
#line 1019 "zlahef_rk.f"
		    zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 1019 "zlahef_rk.f"
		}

/*              Interchange rows KK and KP in first K-1 columns of A */
/*              (column K (or K and K+1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in first KK columns of W. */

#line 1027 "zlahef_rk.f"
		if (k > 1) {
#line 1027 "zlahef_rk.f"
		    i__1 = k - 1;
#line 1027 "zlahef_rk.f"
		    zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 1027 "zlahef_rk.f"
		}
#line 1029 "zlahef_rk.f"
		zswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 1030 "zlahef_rk.f"
	    }

#line 1032 "zlahef_rk.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of W now holds */

/*              W(k) = L(k)*D(k), */

/*              where L(k) is the k-th column of L */

/*              (1) Store subdiag. elements of column L(k) */
/*              and 1-by-1 block D(k) in column k of A. */
/*              (NOTE: Diagonal element L(k,k) is a UNIT element */
/*              and not stored) */
/*                 A(k,k) := D(k,k) = W(k,k) */
/*                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k) */

/*              (NOTE: No need to use for Hermitian matrix */
/*              A( K, K ) = REAL( W( K, K) ) to separately copy diagonal */
/*              element D(k,k) from W (potentially saves only one load)) */
#line 1050 "zlahef_rk.f"
		i__1 = *n - k + 1;
#line 1050 "zlahef_rk.f"
		zcopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 1051 "zlahef_rk.f"
		if (k < *n) {

/*                 (NOTE: No need to check if A(k,k) is NOT ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  case A(k,k) = 0 falls into 2x2 pivot case(3)) */

/*                 Handle division by a small number */

#line 1059 "zlahef_rk.f"
		    i__1 = k + k * a_dim1;
#line 1059 "zlahef_rk.f"
		    t = a[i__1].r;
#line 1060 "zlahef_rk.f"
		    if (abs(t) >= sfmin) {
#line 1061 "zlahef_rk.f"
			r1 = 1. / t;
#line 1062 "zlahef_rk.f"
			i__1 = *n - k;
#line 1062 "zlahef_rk.f"
			zdscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 1063 "zlahef_rk.f"
		    } else {
#line 1064 "zlahef_rk.f"
			i__1 = *n;
#line 1064 "zlahef_rk.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 1065 "zlahef_rk.f"
			    i__2 = ii + k * a_dim1;
#line 1065 "zlahef_rk.f"
			    i__3 = ii + k * a_dim1;
#line 1065 "zlahef_rk.f"
			    z__1.r = a[i__3].r / t, z__1.i = a[i__3].i / t;
#line 1065 "zlahef_rk.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 1066 "zlahef_rk.f"
/* L74: */
#line 1066 "zlahef_rk.f"
			}
#line 1067 "zlahef_rk.f"
		    }

/*                 (2) Conjugate column W(k) */

#line 1071 "zlahef_rk.f"
		    i__1 = *n - k;
#line 1071 "zlahef_rk.f"
		    zlacgv_(&i__1, &w[k + 1 + k * w_dim1], &c__1);

/*                 Store the subdiagonal element of D in array E */

#line 1075 "zlahef_rk.f"
		    i__1 = k;
#line 1075 "zlahef_rk.f"
		    e[i__1].r = 0., e[i__1].i = 0.;

#line 1077 "zlahef_rk.f"
		}

#line 1079 "zlahef_rk.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

/*              (1) Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2 */
/*              block D(k:k+1,k:k+1) in columns k and k+1 of A. */
/*              NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT */
/*              block and not stored. */
/*                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1) */
/*                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) = */
/*                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) ) */

#line 1096 "zlahef_rk.f"
		if (k < *n - 1) {

/*                 Factor out the columns of the inverse of 2-by-2 pivot */
/*                 block D, so that each column contains 1, to reduce the */
/*                 number of FLOPS when we multiply panel */
/*                 ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1). */

/*                 D**(-1) = ( d11 cj(d21) )**(-1) = */
/*                           ( d21    d22 ) */

/*                 = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) = */
/*                                          ( (-d21) (     d11 ) ) */

/*                 = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) * */

/*                   * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) = */
/*                     (     (      -1 )           ( d11/conj(d21) ) ) */

/*                 = 1/(|d21|**2) * 1/(D22*D11-1) * */

/*                   * ( d21*( D11 ) conj(d21)*(  -1 ) ) = */
/*                     (     (  -1 )           ( D22 ) ) */

/*                 = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) = */
/*                                      (     (  -1 )           ( D22 ) ) */

/*                 = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) = */
/*                   (               (  -1 )         ( D22 ) ) */

/*                 Handle division by a small number. (NOTE: order of */
/*                 operations is important) */

/*                 = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) ) */
/*                   (   ((  -1 )          )   (( D22 )     ) ), */

/*                 where D11 = d22/d21, */
/*                       D22 = d11/conj(d21), */
/*                       D21 = d21, */
/*                       T = 1/(D22*D11-1). */

/*                 (NOTE: No need to check for division by ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  (a) d21 != 0 in 2x2 pivot case(4), */
/*                      since |d21| should be larger than |d11| and |d22|; */
/*                  (b) (D22*D11 - 1) != 0, since from (a), */
/*                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.) */

#line 1143 "zlahef_rk.f"
		    i__1 = k + 1 + k * w_dim1;
#line 1143 "zlahef_rk.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 1144 "zlahef_rk.f"
		    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
#line 1144 "zlahef_rk.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 1145 "zlahef_rk.f"
		    d_cnjg(&z__2, &d21);
#line 1145 "zlahef_rk.f"
		    z_div(&z__1, &w[k + k * w_dim1], &z__2);
#line 1145 "zlahef_rk.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 1146 "zlahef_rk.f"
		    z__1.r = d11.r * d22.r - d11.i * d22.i, z__1.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 1146 "zlahef_rk.f"
		    t = 1. / (z__1.r - 1.);

/*                 Update elements in columns A(k) and A(k+1) as */
/*                 dot products of rows of ( W(k) W(k+1) ) and columns */
/*                 of D**(-1) */

#line 1152 "zlahef_rk.f"
		    i__1 = *n;
#line 1152 "zlahef_rk.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 1153 "zlahef_rk.f"
			i__2 = j + k * a_dim1;
#line 1153 "zlahef_rk.f"
			i__3 = j + k * w_dim1;
#line 1153 "zlahef_rk.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 1153 "zlahef_rk.f"
			i__4 = j + (k + 1) * w_dim1;
#line 1153 "zlahef_rk.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 1153 "zlahef_rk.f"
			d_cnjg(&z__5, &d21);
#line 1153 "zlahef_rk.f"
			z_div(&z__2, &z__3, &z__5);
#line 1153 "zlahef_rk.f"
			z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 1153 "zlahef_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 1155 "zlahef_rk.f"
			i__2 = j + (k + 1) * a_dim1;
#line 1155 "zlahef_rk.f"
			i__3 = j + (k + 1) * w_dim1;
#line 1155 "zlahef_rk.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 1155 "zlahef_rk.f"
			i__4 = j + k * w_dim1;
#line 1155 "zlahef_rk.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 1155 "zlahef_rk.f"
			z_div(&z__2, &z__3, &d21);
#line 1155 "zlahef_rk.f"
			z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 1155 "zlahef_rk.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 1157 "zlahef_rk.f"
/* L80: */
#line 1157 "zlahef_rk.f"
		    }
#line 1158 "zlahef_rk.f"
		}

/*              Copy diagonal elements of D(K) to A, */
/*              copy subdiagonal element of D(K) to E(K) and */
/*              ZERO out subdiagonal entry of A */

#line 1164 "zlahef_rk.f"
		i__1 = k + k * a_dim1;
#line 1164 "zlahef_rk.f"
		i__2 = k + k * w_dim1;
#line 1164 "zlahef_rk.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 1165 "zlahef_rk.f"
		i__1 = k + 1 + k * a_dim1;
#line 1165 "zlahef_rk.f"
		a[i__1].r = 0., a[i__1].i = 0.;
#line 1166 "zlahef_rk.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 1166 "zlahef_rk.f"
		i__2 = k + 1 + (k + 1) * w_dim1;
#line 1166 "zlahef_rk.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 1167 "zlahef_rk.f"
		i__1 = k;
#line 1167 "zlahef_rk.f"
		i__2 = k + 1 + k * w_dim1;
#line 1167 "zlahef_rk.f"
		e[i__1].r = w[i__2].r, e[i__1].i = w[i__2].i;
#line 1168 "zlahef_rk.f"
		i__1 = k + 1;
#line 1168 "zlahef_rk.f"
		e[i__1].r = 0., e[i__1].i = 0.;

/*              (2) Conjugate columns W(k) and W(k+1) */

#line 1172 "zlahef_rk.f"
		i__1 = *n - k;
#line 1172 "zlahef_rk.f"
		zlacgv_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 1173 "zlahef_rk.f"
		i__1 = *n - k - 1;
#line 1173 "zlahef_rk.f"
		zlacgv_(&i__1, &w[k + 2 + (k + 1) * w_dim1], &c__1);

#line 1175 "zlahef_rk.f"
	    }

/*           End column K is nonsingular */

#line 1179 "zlahef_rk.f"
	}

/*        Store details of the interchanges in IPIV */

#line 1183 "zlahef_rk.f"
	if (kstep == 1) {
#line 1184 "zlahef_rk.f"
	    ipiv[k] = kp;
#line 1185 "zlahef_rk.f"
	} else {
#line 1186 "zlahef_rk.f"
	    ipiv[k] = -p;
#line 1187 "zlahef_rk.f"
	    ipiv[k + 1] = -kp;
#line 1188 "zlahef_rk.f"
	}

/*        Increase K and return to the start of the main loop */

#line 1192 "zlahef_rk.f"
	k += kstep;
#line 1193 "zlahef_rk.f"
	goto L70;

#line 1195 "zlahef_rk.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**H = A22 - L21*W**H */

/*        computing blocks of NB columns at a time (note that conjg(W) is */
/*        actually stored) */

#line 1204 "zlahef_rk.f"
	i__1 = *n;
#line 1204 "zlahef_rk.f"
	i__2 = *nb;
#line 1204 "zlahef_rk.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 1205 "zlahef_rk.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 1205 "zlahef_rk.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 1209 "zlahef_rk.f"
	    i__3 = j + jb - 1;
#line 1209 "zlahef_rk.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 1210 "zlahef_rk.f"
		i__4 = jj + jj * a_dim1;
#line 1210 "zlahef_rk.f"
		i__5 = jj + jj * a_dim1;
#line 1210 "zlahef_rk.f"
		d__1 = a[i__5].r;
#line 1210 "zlahef_rk.f"
		a[i__4].r = d__1, a[i__4].i = 0.;
#line 1211 "zlahef_rk.f"
		i__4 = j + jb - jj;
#line 1211 "zlahef_rk.f"
		i__5 = k - 1;
#line 1211 "zlahef_rk.f"
		z__1.r = -1., z__1.i = -0.;
#line 1211 "zlahef_rk.f"
		zgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1]
			, &c__1, (ftnlen)12);
#line 1214 "zlahef_rk.f"
		i__4 = jj + jj * a_dim1;
#line 1214 "zlahef_rk.f"
		i__5 = jj + jj * a_dim1;
#line 1214 "zlahef_rk.f"
		d__1 = a[i__5].r;
#line 1214 "zlahef_rk.f"
		a[i__4].r = d__1, a[i__4].i = 0.;
#line 1215 "zlahef_rk.f"
/* L100: */
#line 1215 "zlahef_rk.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 1219 "zlahef_rk.f"
	    if (j + jb <= *n) {
#line 1219 "zlahef_rk.f"
		i__3 = *n - j - jb + 1;
#line 1219 "zlahef_rk.f"
		i__4 = k - 1;
#line 1219 "zlahef_rk.f"
		z__1.r = -1., z__1.i = -0.;
#line 1219 "zlahef_rk.f"
		zgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, 
			&a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 1219 "zlahef_rk.f"
	    }
#line 1223 "zlahef_rk.f"
/* L110: */
#line 1223 "zlahef_rk.f"
	}

/*        Set KB to the number of columns factorized */

#line 1227 "zlahef_rk.f"
	*kb = k - 1;

#line 1229 "zlahef_rk.f"
    }
#line 1230 "zlahef_rk.f"
    return 0;

/*     End of ZLAHEF_RK */

} /* zlahef_rk__ */

