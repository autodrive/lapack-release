#line 1 "zsytrf_rk.f"
/* zsytrf_rk.f -- translated by f2c (version 20100827).
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

#line 1 "zsytrf_rk.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b ZSYTRF_RK computes the factorization of a complex symmetric indefinite matrix using the bounded
 Bunch-Kaufman (rook) diagonal pivoting method (BLAS3 blocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRF_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrf_
rk.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrf_
rk.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrf_
rk.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, */
/*                             INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), E ( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > ZSYTRF_RK computes the factorization of a complex symmetric matrix A */
/* > using the bounded Bunch-Kaufman (rook) diagonal pivoting method: */
/* > */
/* >    A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is symmetric and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension ( MAX(1,LWORK) ). */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of WORK.  LWORK >=1.  For best performance */
/* >          LWORK >= N*NB, where NB is the block size returned */
/* >          by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; */
/* >          the routine only calculates the optimal size of the WORK */
/* >          array, returns this value as the first entry of the WORK */
/* >          array, and no error message related to LWORK is issued */
/* >          by XERBLA. */
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
/* > TODO: put correct description */
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
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zsytrf_rk__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *e, integer *ipiv, doublecomplex *work, 
	integer *lwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int zsytf2_rk__(char *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, integer *, ftnlen), 
	    zlasyf_rk__(char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer kb, nb, ip, iws;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

#line 298 "zsytrf_rk.f"
    /* Parameter adjustments */
#line 298 "zsytrf_rk.f"
    a_dim1 = *lda;
#line 298 "zsytrf_rk.f"
    a_offset = 1 + a_dim1;
#line 298 "zsytrf_rk.f"
    a -= a_offset;
#line 298 "zsytrf_rk.f"
    --e;
#line 298 "zsytrf_rk.f"
    --ipiv;
#line 298 "zsytrf_rk.f"
    --work;
#line 298 "zsytrf_rk.f"

#line 298 "zsytrf_rk.f"
    /* Function Body */
#line 298 "zsytrf_rk.f"
    *info = 0;
#line 299 "zsytrf_rk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 300 "zsytrf_rk.f"
    lquery = *lwork == -1;
#line 301 "zsytrf_rk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 302 "zsytrf_rk.f"
	*info = -1;
#line 303 "zsytrf_rk.f"
    } else if (*n < 0) {
#line 304 "zsytrf_rk.f"
	*info = -2;
#line 305 "zsytrf_rk.f"
    } else if (*lda < max(1,*n)) {
#line 306 "zsytrf_rk.f"
	*info = -4;
#line 307 "zsytrf_rk.f"
    } else if (*lwork < 1 && ! lquery) {
#line 308 "zsytrf_rk.f"
	*info = -8;
#line 309 "zsytrf_rk.f"
    }

#line 311 "zsytrf_rk.f"
    if (*info == 0) {

/*        Determine the block size */

#line 315 "zsytrf_rk.f"
	nb = ilaenv_(&c__1, "ZSYTRF_RK", uplo, n, &c_n1, &c_n1, &c_n1, (
		ftnlen)9, (ftnlen)1);
#line 316 "zsytrf_rk.f"
	lwkopt = *n * nb;
#line 317 "zsytrf_rk.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 318 "zsytrf_rk.f"
    }

#line 320 "zsytrf_rk.f"
    if (*info != 0) {
#line 321 "zsytrf_rk.f"
	i__1 = -(*info);
#line 321 "zsytrf_rk.f"
	xerbla_("ZSYTRF_RK", &i__1, (ftnlen)9);
#line 322 "zsytrf_rk.f"
	return 0;
#line 323 "zsytrf_rk.f"
    } else if (lquery) {
#line 324 "zsytrf_rk.f"
	return 0;
#line 325 "zsytrf_rk.f"
    }

#line 327 "zsytrf_rk.f"
    nbmin = 2;
#line 328 "zsytrf_rk.f"
    ldwork = *n;
#line 329 "zsytrf_rk.f"
    if (nb > 1 && nb < *n) {
#line 330 "zsytrf_rk.f"
	iws = ldwork * nb;
#line 331 "zsytrf_rk.f"
	if (*lwork < iws) {
/* Computing MAX */
#line 332 "zsytrf_rk.f"
	    i__1 = *lwork / ldwork;
#line 332 "zsytrf_rk.f"
	    nb = max(i__1,1);
/* Computing MAX */
#line 333 "zsytrf_rk.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "ZSYTRF_RK", uplo, n, &c_n1, &
		    c_n1, &c_n1, (ftnlen)9, (ftnlen)1);
#line 333 "zsytrf_rk.f"
	    nbmin = max(i__1,i__2);
#line 335 "zsytrf_rk.f"
	}
#line 336 "zsytrf_rk.f"
    } else {
#line 337 "zsytrf_rk.f"
	iws = 1;
#line 338 "zsytrf_rk.f"
    }
#line 339 "zsytrf_rk.f"
    if (nb < nbmin) {
#line 339 "zsytrf_rk.f"
	nb = *n;
#line 339 "zsytrf_rk.f"
    }

#line 342 "zsytrf_rk.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        KB, where KB is the number of columns factorized by ZLASYF_RK; */
/*        KB is either NB or NB-1, or K for the last block */

#line 350 "zsytrf_rk.f"
	k = *n;
#line 351 "zsytrf_rk.f"
L10:

/*        If K < 1, exit from loop */

#line 355 "zsytrf_rk.f"
	if (k < 1) {
#line 355 "zsytrf_rk.f"
	    goto L15;
#line 355 "zsytrf_rk.f"
	}

#line 358 "zsytrf_rk.f"
	if (k > nb) {

/*           Factorize columns k-kb+1:k of A and use blocked code to */
/*           update columns 1:k-kb */

#line 363 "zsytrf_rk.f"
	    zlasyf_rk__(uplo, &k, &nb, &kb, &a[a_offset], lda, &e[1], &ipiv[1]
		    , &work[1], &ldwork, &iinfo, (ftnlen)1);
#line 365 "zsytrf_rk.f"
	} else {

/*           Use unblocked code to factorize columns 1:k of A */

#line 369 "zsytrf_rk.f"
	    zsytf2_rk__(uplo, &k, &a[a_offset], lda, &e[1], &ipiv[1], &iinfo, 
		    (ftnlen)1);
#line 370 "zsytrf_rk.f"
	    kb = k;
#line 371 "zsytrf_rk.f"
	}

/*        Set INFO on the first occurrence of a zero pivot */

#line 375 "zsytrf_rk.f"
	if (*info == 0 && iinfo > 0) {
#line 375 "zsytrf_rk.f"
	    *info = iinfo;
#line 375 "zsytrf_rk.f"
	}

/*        No need to adjust IPIV */


/*        Apply permutations to the leading panel 1:k-1 */

/*        Read IPIV from the last block factored, i.e. */
/*        indices  k-kb+1:k and apply row permutations to the */
/*        last k+1 colunms k+1:N after that block */
/*        (We can do the simple loop over IPIV with decrement -1, */
/*        since the ABS value of IPIV( I ) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 390 "zsytrf_rk.f"
	if (k < *n) {
#line 391 "zsytrf_rk.f"
	    i__1 = k - kb + 1;
#line 391 "zsytrf_rk.f"
	    for (i__ = k; i__ >= i__1; --i__) {
#line 392 "zsytrf_rk.f"
		ip = (i__2 = ipiv[i__], abs(i__2));
#line 393 "zsytrf_rk.f"
		if (ip != i__) {
#line 394 "zsytrf_rk.f"
		    i__2 = *n - k;
#line 394 "zsytrf_rk.f"
		    zswap_(&i__2, &a[i__ + (k + 1) * a_dim1], lda, &a[ip + (k 
			    + 1) * a_dim1], lda);
#line 396 "zsytrf_rk.f"
		}
#line 397 "zsytrf_rk.f"
	    }
#line 398 "zsytrf_rk.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 402 "zsytrf_rk.f"
	k -= kb;
#line 403 "zsytrf_rk.f"
	goto L10;

/*        This label is the exit from main loop over K decreasing */
/*        from N to 1 in steps of KB */

#line 408 "zsytrf_rk.f"
L15:

#line 410 "zsytrf_rk.f"
	;
#line 410 "zsytrf_rk.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        KB, where KB is the number of columns factorized by ZLASYF_RK; */
/*        KB is either NB or NB-1, or N-K+1 for the last block */

#line 418 "zsytrf_rk.f"
	k = 1;
#line 419 "zsytrf_rk.f"
L20:

/*        If K > N, exit from loop */

#line 423 "zsytrf_rk.f"
	if (k > *n) {
#line 423 "zsytrf_rk.f"
	    goto L35;
#line 423 "zsytrf_rk.f"
	}

#line 426 "zsytrf_rk.f"
	if (k <= *n - nb) {

/*           Factorize columns k:k+kb-1 of A and use blocked code to */
/*           update columns k+kb:n */

#line 431 "zsytrf_rk.f"
	    i__1 = *n - k + 1;
#line 431 "zsytrf_rk.f"
	    zlasyf_rk__(uplo, &i__1, &nb, &kb, &a[k + k * a_dim1], lda, &e[k],
		     &ipiv[k], &work[1], &ldwork, &iinfo, (ftnlen)1);
#line 435 "zsytrf_rk.f"
	} else {

/*           Use unblocked code to factorize columns k:n of A */

#line 439 "zsytrf_rk.f"
	    i__1 = *n - k + 1;
#line 439 "zsytrf_rk.f"
	    zsytf2_rk__(uplo, &i__1, &a[k + k * a_dim1], lda, &e[k], &ipiv[k],
		     &iinfo, (ftnlen)1);
#line 441 "zsytrf_rk.f"
	    kb = *n - k + 1;

#line 443 "zsytrf_rk.f"
	}

/*        Set INFO on the first occurrence of a zero pivot */

#line 447 "zsytrf_rk.f"
	if (*info == 0 && iinfo > 0) {
#line 447 "zsytrf_rk.f"
	    *info = iinfo + k - 1;
#line 447 "zsytrf_rk.f"
	}

/*        Adjust IPIV */

#line 452 "zsytrf_rk.f"
	i__1 = k + kb - 1;
#line 452 "zsytrf_rk.f"
	for (i__ = k; i__ <= i__1; ++i__) {
#line 453 "zsytrf_rk.f"
	    if (ipiv[i__] > 0) {
#line 454 "zsytrf_rk.f"
		ipiv[i__] = ipiv[i__] + k - 1;
#line 455 "zsytrf_rk.f"
	    } else {
#line 456 "zsytrf_rk.f"
		ipiv[i__] = ipiv[i__] - k + 1;
#line 457 "zsytrf_rk.f"
	    }
#line 458 "zsytrf_rk.f"
	}

/*        Apply permutations to the leading panel 1:k-1 */

/*        Read IPIV from the last block factored, i.e. */
/*        indices  k:k+kb-1 and apply row permutations to the */
/*        first k-1 colunms 1:k-1 before that block */
/*        (We can do the simple loop over IPIV with increment 1, */
/*        since the ABS value of IPIV( I ) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 469 "zsytrf_rk.f"
	if (k > 1) {
#line 470 "zsytrf_rk.f"
	    i__1 = k + kb - 1;
#line 470 "zsytrf_rk.f"
	    for (i__ = k; i__ <= i__1; ++i__) {
#line 471 "zsytrf_rk.f"
		ip = (i__2 = ipiv[i__], abs(i__2));
#line 472 "zsytrf_rk.f"
		if (ip != i__) {
#line 473 "zsytrf_rk.f"
		    i__2 = k - 1;
#line 473 "zsytrf_rk.f"
		    zswap_(&i__2, &a[i__ + a_dim1], lda, &a[ip + a_dim1], lda)
			    ;
#line 475 "zsytrf_rk.f"
		}
#line 476 "zsytrf_rk.f"
	    }
#line 477 "zsytrf_rk.f"
	}

/*        Increase K and return to the start of the main loop */

#line 481 "zsytrf_rk.f"
	k += kb;
#line 482 "zsytrf_rk.f"
	goto L20;

/*        This label is the exit from main loop over K increasing */
/*        from 1 to N in steps of KB */

#line 487 "zsytrf_rk.f"
L35:

/*     End Lower */

#line 491 "zsytrf_rk.f"
	;
#line 491 "zsytrf_rk.f"
    }

#line 493 "zsytrf_rk.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 494 "zsytrf_rk.f"
    return 0;

/*     End of ZSYTRF_RK */

} /* zsytrf_rk__ */

