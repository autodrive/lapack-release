#line 1 "slasyf_rook.f"
/* slasyf_rook.f -- translated by f2c (version 20100827).
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

#line 1 "slasyf_rook.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = -1.;
static doublereal c_b10 = 1.;

/* > \brief \b SLASYF_ROOK computes a partial factorization of a real symmetric matrix using the bounded Bunch
-Kaufman ("rook") diagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASYF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASYF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KB, LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASYF_ROOK computes a partial factorization of a real symmetric */
/* > matrix A using the bounded Bunch-Kaufman ("rook") diagonal */
/* > pivoting method. The partial factorization has the form: */
/* > */
/* > A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or: */
/* >       ( 0  U22 ) (  0   D  ) ( U12**T U22**T ) */
/* > */
/* > A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L' */
/* >       ( L21  I ) (  0  A22 ) (  0       I    ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > */
/* > SLASYF_ROOK is an auxiliary routine called by SSYTRF_ROOK. It uses */
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
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          n-by-n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n-by-n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit, A contains details of the partial factorization. */
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
/* >             Only the last KB elements of IPIV are set. */
/* > */
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >             interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and */
/* >             columns k and -IPIV(k) were interchanged and rows and */
/* >             columns k-1 and -IPIV(k-1) were inerchaged, */
/* >             D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
/* > */
/* >          If UPLO = 'L': */
/* >             Only the first KB elements of IPIV are set. */
/* > */
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) */
/* >             were interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and */
/* >             columns k and -IPIV(k) were interchanged and rows and */
/* >             columns k+1 and -IPIV(k+1) were inerchaged, */
/* >             D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
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
/* >          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization */
/* >               has been completed, but the block diagonal matrix D is */
/* >               exactly singular. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2013 */

/* > \ingroup realSYcomputational */

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
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int slasyf_rook__(char *uplo, integer *n, integer *nb, 
	integer *kb, doublereal *a, integer *lda, integer *ipiv, doublereal *
	w, integer *ldw, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, p;
    static doublereal t, r1, d11, d12, d21, d22;
    static integer jb, ii, jj, kk, kp, kw, jp1, jp2, kkw;
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

#line 231 "slasyf_rook.f"
    /* Parameter adjustments */
#line 231 "slasyf_rook.f"
    a_dim1 = *lda;
#line 231 "slasyf_rook.f"
    a_offset = 1 + a_dim1;
#line 231 "slasyf_rook.f"
    a -= a_offset;
#line 231 "slasyf_rook.f"
    --ipiv;
#line 231 "slasyf_rook.f"
    w_dim1 = *ldw;
#line 231 "slasyf_rook.f"
    w_offset = 1 + w_dim1;
#line 231 "slasyf_rook.f"
    w -= w_offset;
#line 231 "slasyf_rook.f"

#line 231 "slasyf_rook.f"
    /* Function Body */
#line 231 "slasyf_rook.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 235 "slasyf_rook.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 239 "slasyf_rook.f"
    sfmin = slamch_("S", (ftnlen)1);

#line 241 "slasyf_rook.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 */

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

#line 249 "slasyf_rook.f"
	k = *n;
#line 250 "slasyf_rook.f"
L10:

/*        KW is the column of W which corresponds to column K of A */

#line 254 "slasyf_rook.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 258 "slasyf_rook.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 258 "slasyf_rook.f"
	    goto L30;
#line 258 "slasyf_rook.f"
	}

#line 261 "slasyf_rook.f"
	kstep = 1;
#line 262 "slasyf_rook.f"
	p = k;

/*        Copy column K of A to column KW of W and update it */

#line 266 "slasyf_rook.f"
	scopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
#line 267 "slasyf_rook.f"
	if (k < *n) {
#line 267 "slasyf_rook.f"
	    i__1 = *n - k;
#line 267 "slasyf_rook.f"
	    sgemv_("No transpose", &k, &i__1, &c_b9, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b10, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 267 "slasyf_rook.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 274 "slasyf_rook.f"
	absakk = (d__1 = w[k + kw * w_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 280 "slasyf_rook.f"
	if (k > 1) {
#line 281 "slasyf_rook.f"
	    i__1 = k - 1;
#line 281 "slasyf_rook.f"
	    imax = isamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 282 "slasyf_rook.f"
	    colmax = (d__1 = w[imax + kw * w_dim1], abs(d__1));
#line 283 "slasyf_rook.f"
	} else {
#line 284 "slasyf_rook.f"
	    colmax = 0.;
#line 285 "slasyf_rook.f"
	}

#line 287 "slasyf_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 291 "slasyf_rook.f"
	    if (*info == 0) {
#line 291 "slasyf_rook.f"
		*info = k;
#line 291 "slasyf_rook.f"
	    }
#line 293 "slasyf_rook.f"
	    kp = k;
#line 294 "slasyf_rook.f"
	    scopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
#line 295 "slasyf_rook.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 304 "slasyf_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 308 "slasyf_rook.f"
		kp = k;

#line 310 "slasyf_rook.f"
	    } else {

#line 312 "slasyf_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 316 "slasyf_rook.f"
L12:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column KW-1 of W and update it */

#line 323 "slasyf_rook.f"
		scopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
#line 324 "slasyf_rook.f"
		i__1 = k - imax;
#line 324 "slasyf_rook.f"
		scopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);

#line 327 "slasyf_rook.f"
		if (k < *n) {
#line 327 "slasyf_rook.f"
		    i__1 = *n - k;
#line 327 "slasyf_rook.f"
		    sgemv_("No transpose", &k, &i__1, &c_b9, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b10, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 327 "slasyf_rook.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 336 "slasyf_rook.f"
		if (imax != k) {
#line 337 "slasyf_rook.f"
		    i__1 = k - imax;
#line 337 "slasyf_rook.f"
		    jmax = imax + isamax_(&i__1, &w[imax + 1 + (kw - 1) * 
			    w_dim1], &c__1);
#line 339 "slasyf_rook.f"
		    rowmax = (d__1 = w[jmax + (kw - 1) * w_dim1], abs(d__1));
#line 340 "slasyf_rook.f"
		} else {
#line 341 "slasyf_rook.f"
		    rowmax = 0.;
#line 342 "slasyf_rook.f"
		}

#line 344 "slasyf_rook.f"
		if (imax > 1) {
#line 345 "slasyf_rook.f"
		    i__1 = imax - 1;
#line 345 "slasyf_rook.f"
		    itemp = isamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
#line 346 "slasyf_rook.f"
		    stemp = (d__1 = w[itemp + (kw - 1) * w_dim1], abs(d__1));
#line 347 "slasyf_rook.f"
		    if (stemp > rowmax) {
#line 348 "slasyf_rook.f"
			rowmax = stemp;
#line 349 "slasyf_rook.f"
			jmax = itemp;
#line 350 "slasyf_rook.f"
		    }
#line 351 "slasyf_rook.f"
		}

/*                 Equivalent to testing for */
/*                 ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 357 "slasyf_rook.f"
		if (! ((d__1 = w[imax + (kw - 1) * w_dim1], abs(d__1)) < 
			alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 363 "slasyf_rook.f"
		    kp = imax;

/*                    copy column KW-1 of W to column KW of W */

#line 367 "slasyf_rook.f"
		    scopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 369 "slasyf_rook.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 374 "slasyf_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 380 "slasyf_rook.f"
		    kp = imax;
#line 381 "slasyf_rook.f"
		    kstep = 2;
#line 382 "slasyf_rook.f"
		    done = TRUE_;
#line 383 "slasyf_rook.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 387 "slasyf_rook.f"
		    p = imax;
#line 388 "slasyf_rook.f"
		    colmax = rowmax;
#line 389 "slasyf_rook.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 393 "slasyf_rook.f"
		    scopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 395 "slasyf_rook.f"
		}

/*                 End pivot search loop body */

#line 399 "slasyf_rook.f"
		if (! done) {
#line 399 "slasyf_rook.f"
		    goto L12;
#line 399 "slasyf_rook.f"
		}

#line 401 "slasyf_rook.f"
	    }

/*           ============================================================ */

#line 405 "slasyf_rook.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 409 "slasyf_rook.f"
	    kkw = *nb + kk - *n;

#line 411 "slasyf_rook.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

#line 415 "slasyf_rook.f"
		i__1 = k - p;
#line 415 "slasyf_rook.f"
		scopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * 
			a_dim1], lda);
#line 416 "slasyf_rook.f"
		scopy_(&p, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &
			c__1);

/*              Interchange rows K and P in last N-K+1 columns of A */
/*              and last N-K+2 columns of W */

#line 421 "slasyf_rook.f"
		i__1 = *n - k + 1;
#line 421 "slasyf_rook.f"
		sswap_(&i__1, &a[k + k * a_dim1], lda, &a[p + k * a_dim1], 
			lda);
#line 422 "slasyf_rook.f"
		i__1 = *n - kk + 1;
#line 422 "slasyf_rook.f"
		sswap_(&i__1, &w[k + kkw * w_dim1], ldw, &w[p + kkw * w_dim1],
			 ldw);
#line 423 "slasyf_rook.f"
	    }

/*           Updated column KP is already stored in column KKW of W */

#line 427 "slasyf_rook.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

#line 431 "slasyf_rook.f"
		a[kp + k * a_dim1] = a[kk + k * a_dim1];
#line 432 "slasyf_rook.f"
		i__1 = k - 1 - kp;
#line 432 "slasyf_rook.f"
		scopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 434 "slasyf_rook.f"
		scopy_(&kp, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
			c__1);

/*              Interchange rows KK and KP in last N-KK+1 columns */
/*              of A and W */

#line 439 "slasyf_rook.f"
		i__1 = *n - kk + 1;
#line 439 "slasyf_rook.f"
		sswap_(&i__1, &a[kk + kk * a_dim1], lda, &a[kp + kk * a_dim1],
			 lda);
#line 440 "slasyf_rook.f"
		i__1 = *n - kk + 1;
#line 440 "slasyf_rook.f"
		sswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 442 "slasyf_rook.f"
	    }

#line 444 "slasyf_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column KW of W now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Store U(k) in column k of A */

#line 454 "slasyf_rook.f"
		scopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 455 "slasyf_rook.f"
		if (k > 1) {
#line 456 "slasyf_rook.f"
		    if ((d__1 = a[k + k * a_dim1], abs(d__1)) >= sfmin) {
#line 457 "slasyf_rook.f"
			r1 = 1. / a[k + k * a_dim1];
#line 458 "slasyf_rook.f"
			i__1 = k - 1;
#line 458 "slasyf_rook.f"
			sscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 459 "slasyf_rook.f"
		    } else if (a[k + k * a_dim1] != 0.) {
#line 460 "slasyf_rook.f"
			i__1 = k - 1;
#line 460 "slasyf_rook.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 461 "slasyf_rook.f"
			    a[ii + k * a_dim1] /= a[k + k * a_dim1];
#line 462 "slasyf_rook.f"
/* L14: */
#line 462 "slasyf_rook.f"
			}
#line 463 "slasyf_rook.f"
		    }
#line 464 "slasyf_rook.f"
		}

#line 466 "slasyf_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns KW and KW-1 of W now */
/*              hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

#line 476 "slasyf_rook.f"
		if (k > 2) {

/*                 Store U(k) and U(k-1) in columns k and k-1 of A */

#line 480 "slasyf_rook.f"
		    d12 = w[k - 1 + kw * w_dim1];
#line 481 "slasyf_rook.f"
		    d11 = w[k + kw * w_dim1] / d12;
#line 482 "slasyf_rook.f"
		    d22 = w[k - 1 + (kw - 1) * w_dim1] / d12;
#line 483 "slasyf_rook.f"
		    t = 1. / (d11 * d22 - 1.);
#line 484 "slasyf_rook.f"
		    i__1 = k - 2;
#line 484 "slasyf_rook.f"
		    for (j = 1; j <= i__1; ++j) {
#line 485 "slasyf_rook.f"
			a[j + (k - 1) * a_dim1] = t * ((d11 * w[j + (kw - 1) *
				 w_dim1] - w[j + kw * w_dim1]) / d12);
#line 487 "slasyf_rook.f"
			a[j + k * a_dim1] = t * ((d22 * w[j + kw * w_dim1] - 
				w[j + (kw - 1) * w_dim1]) / d12);
#line 489 "slasyf_rook.f"
/* L20: */
#line 489 "slasyf_rook.f"
		    }
#line 490 "slasyf_rook.f"
		}

/*              Copy D(k) to A */

#line 494 "slasyf_rook.f"
		a[k - 1 + (k - 1) * a_dim1] = w[k - 1 + (kw - 1) * w_dim1];
#line 495 "slasyf_rook.f"
		a[k - 1 + k * a_dim1] = w[k - 1 + kw * w_dim1];
#line 496 "slasyf_rook.f"
		a[k + k * a_dim1] = w[k + kw * w_dim1];
#line 497 "slasyf_rook.f"
	    }
#line 498 "slasyf_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 502 "slasyf_rook.f"
	if (kstep == 1) {
#line 503 "slasyf_rook.f"
	    ipiv[k] = kp;
#line 504 "slasyf_rook.f"
	} else {
#line 505 "slasyf_rook.f"
	    ipiv[k] = -p;
#line 506 "slasyf_rook.f"
	    ipiv[k - 1] = -kp;
#line 507 "slasyf_rook.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 511 "slasyf_rook.f"
	k -= kstep;
#line 512 "slasyf_rook.f"
	goto L10;

#line 514 "slasyf_rook.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T */

/*        computing blocks of NB columns at a time */

#line 522 "slasyf_rook.f"
	i__1 = -(*nb);
#line 522 "slasyf_rook.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 523 "slasyf_rook.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 523 "slasyf_rook.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 527 "slasyf_rook.f"
	    i__2 = j + jb - 1;
#line 527 "slasyf_rook.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 528 "slasyf_rook.f"
		i__3 = jj - j + 1;
#line 528 "slasyf_rook.f"
		i__4 = *n - k;
#line 528 "slasyf_rook.f"
		sgemv_("No transpose", &i__3, &i__4, &c_b9, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b10,
			 &a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 531 "slasyf_rook.f"
/* L40: */
#line 531 "slasyf_rook.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 535 "slasyf_rook.f"
	    if (j >= 2) {
#line 535 "slasyf_rook.f"
		i__2 = j - 1;
#line 535 "slasyf_rook.f"
		i__3 = *n - k;
#line 535 "slasyf_rook.f"
		sgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &c_b9, 
			&a[(k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * 
			w_dim1], ldw, &c_b10, &a[j * a_dim1 + 1], lda, (
			ftnlen)12, (ftnlen)9);
#line 535 "slasyf_rook.f"
	    }
#line 539 "slasyf_rook.f"
/* L50: */
#line 539 "slasyf_rook.f"
	}

/*        Put U12 in standard form by partially undoing the interchanges */
/*        in columns k+1:n */

#line 544 "slasyf_rook.f"
	j = k + 1;
#line 545 "slasyf_rook.f"
L60:

#line 547 "slasyf_rook.f"
	kstep = 1;
#line 548 "slasyf_rook.f"
	jp1 = 1;
#line 549 "slasyf_rook.f"
	jj = j;
#line 550 "slasyf_rook.f"
	jp2 = ipiv[j];
#line 551 "slasyf_rook.f"
	if (jp2 < 0) {
#line 552 "slasyf_rook.f"
	    jp2 = -jp2;
#line 553 "slasyf_rook.f"
	    ++j;
#line 554 "slasyf_rook.f"
	    jp1 = -ipiv[j];
#line 555 "slasyf_rook.f"
	    kstep = 2;
#line 556 "slasyf_rook.f"
	}

#line 558 "slasyf_rook.f"
	++j;
#line 559 "slasyf_rook.f"
	if (jp2 != jj && j <= *n) {
#line 559 "slasyf_rook.f"
	    i__1 = *n - j + 1;
#line 559 "slasyf_rook.f"
	    sswap_(&i__1, &a[jp2 + j * a_dim1], lda, &a[jj + j * a_dim1], lda)
		    ;
#line 559 "slasyf_rook.f"
	}
#line 561 "slasyf_rook.f"
	jj = j - 1;
#line 562 "slasyf_rook.f"
	if (jp1 != jj && kstep == 2) {
#line 562 "slasyf_rook.f"
	    i__1 = *n - j + 1;
#line 562 "slasyf_rook.f"
	    sswap_(&i__1, &a[jp1 + j * a_dim1], lda, &a[jj + j * a_dim1], lda)
		    ;
#line 562 "slasyf_rook.f"
	}
#line 564 "slasyf_rook.f"
	if (j <= *n) {
#line 564 "slasyf_rook.f"
	    goto L60;
#line 564 "slasyf_rook.f"
	}

/*        Set KB to the number of columns factorized */

#line 569 "slasyf_rook.f"
	*kb = *n - k;

#line 571 "slasyf_rook.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 */

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 579 "slasyf_rook.f"
	k = 1;
#line 580 "slasyf_rook.f"
L70:

/*        Exit from loop */

#line 584 "slasyf_rook.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 584 "slasyf_rook.f"
	    goto L90;
#line 584 "slasyf_rook.f"
	}

#line 587 "slasyf_rook.f"
	kstep = 1;
#line 588 "slasyf_rook.f"
	p = k;

/*        Copy column K of A to column K of W and update it */

#line 592 "slasyf_rook.f"
	i__1 = *n - k + 1;
#line 592 "slasyf_rook.f"
	scopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
#line 593 "slasyf_rook.f"
	if (k > 1) {
#line 593 "slasyf_rook.f"
	    i__1 = *n - k + 1;
#line 593 "slasyf_rook.f"
	    i__2 = k - 1;
#line 593 "slasyf_rook.f"
	    sgemv_("No transpose", &i__1, &i__2, &c_b9, &a[k + a_dim1], lda, &
		    w[k + w_dim1], ldw, &c_b10, &w[k + k * w_dim1], &c__1, (
		    ftnlen)12);
#line 593 "slasyf_rook.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 600 "slasyf_rook.f"
	absakk = (d__1 = w[k + k * w_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 606 "slasyf_rook.f"
	if (k < *n) {
#line 607 "slasyf_rook.f"
	    i__1 = *n - k;
#line 607 "slasyf_rook.f"
	    imax = k + isamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 608 "slasyf_rook.f"
	    colmax = (d__1 = w[imax + k * w_dim1], abs(d__1));
#line 609 "slasyf_rook.f"
	} else {
#line 610 "slasyf_rook.f"
	    colmax = 0.;
#line 611 "slasyf_rook.f"
	}

#line 613 "slasyf_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 617 "slasyf_rook.f"
	    if (*info == 0) {
#line 617 "slasyf_rook.f"
		*info = k;
#line 617 "slasyf_rook.f"
	    }
#line 619 "slasyf_rook.f"
	    kp = k;
#line 620 "slasyf_rook.f"
	    i__1 = *n - k + 1;
#line 620 "slasyf_rook.f"
	    scopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
		    c__1);
#line 621 "slasyf_rook.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 630 "slasyf_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 634 "slasyf_rook.f"
		kp = k;

#line 636 "slasyf_rook.f"
	    } else {

#line 638 "slasyf_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 642 "slasyf_rook.f"
L72:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column K+1 of W and update it */

#line 649 "slasyf_rook.f"
		i__1 = imax - k;
#line 649 "slasyf_rook.f"
		scopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 650 "slasyf_rook.f"
		i__1 = *n - imax + 1;
#line 650 "slasyf_rook.f"
		scopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 
			1) * w_dim1], &c__1);
#line 652 "slasyf_rook.f"
		if (k > 1) {
#line 652 "slasyf_rook.f"
		    i__1 = *n - k + 1;
#line 652 "slasyf_rook.f"
		    i__2 = k - 1;
#line 652 "slasyf_rook.f"
		    sgemv_("No transpose", &i__1, &i__2, &c_b9, &a[k + a_dim1]
			    , lda, &w[imax + w_dim1], ldw, &c_b10, &w[k + (k 
			    + 1) * w_dim1], &c__1, (ftnlen)12);
#line 652 "slasyf_rook.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 661 "slasyf_rook.f"
		if (imax != k) {
#line 662 "slasyf_rook.f"
		    i__1 = imax - k;
#line 662 "slasyf_rook.f"
		    jmax = k - 1 + isamax_(&i__1, &w[k + (k + 1) * w_dim1], &
			    c__1);
#line 663 "slasyf_rook.f"
		    rowmax = (d__1 = w[jmax + (k + 1) * w_dim1], abs(d__1));
#line 664 "slasyf_rook.f"
		} else {
#line 665 "slasyf_rook.f"
		    rowmax = 0.;
#line 666 "slasyf_rook.f"
		}

#line 668 "slasyf_rook.f"
		if (imax < *n) {
#line 669 "slasyf_rook.f"
		    i__1 = *n - imax;
#line 669 "slasyf_rook.f"
		    itemp = imax + isamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
#line 670 "slasyf_rook.f"
		    stemp = (d__1 = w[itemp + (k + 1) * w_dim1], abs(d__1));
#line 671 "slasyf_rook.f"
		    if (stemp > rowmax) {
#line 672 "slasyf_rook.f"
			rowmax = stemp;
#line 673 "slasyf_rook.f"
			jmax = itemp;
#line 674 "slasyf_rook.f"
		    }
#line 675 "slasyf_rook.f"
		}

/*                 Equivalent to testing for */
/*                 ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 681 "slasyf_rook.f"
		if (! ((d__1 = w[imax + (k + 1) * w_dim1], abs(d__1)) < alpha 
			* rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 687 "slasyf_rook.f"
		    kp = imax;

/*                    copy column K+1 of W to column K of W */

#line 691 "slasyf_rook.f"
		    i__1 = *n - k + 1;
#line 691 "slasyf_rook.f"
		    scopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 693 "slasyf_rook.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 698 "slasyf_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 704 "slasyf_rook.f"
		    kp = imax;
#line 705 "slasyf_rook.f"
		    kstep = 2;
#line 706 "slasyf_rook.f"
		    done = TRUE_;
#line 707 "slasyf_rook.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 711 "slasyf_rook.f"
		    p = imax;
#line 712 "slasyf_rook.f"
		    colmax = rowmax;
#line 713 "slasyf_rook.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 717 "slasyf_rook.f"
		    i__1 = *n - k + 1;
#line 717 "slasyf_rook.f"
		    scopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 719 "slasyf_rook.f"
		}

/*                 End pivot search loop body */

#line 723 "slasyf_rook.f"
		if (! done) {
#line 723 "slasyf_rook.f"
		    goto L72;
#line 723 "slasyf_rook.f"
		}

#line 725 "slasyf_rook.f"
	    }

/*           ============================================================ */

#line 729 "slasyf_rook.f"
	    kk = k + kstep - 1;

#line 731 "slasyf_rook.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

#line 735 "slasyf_rook.f"
		i__1 = p - k;
#line 735 "slasyf_rook.f"
		scopy_(&i__1, &a[k + k * a_dim1], &c__1, &a[p + k * a_dim1], 
			lda);
#line 736 "slasyf_rook.f"
		i__1 = *n - p + 1;
#line 736 "slasyf_rook.f"
		scopy_(&i__1, &a[p + k * a_dim1], &c__1, &a[p + p * a_dim1], &
			c__1);

/*              Interchange rows K and P in first K columns of A */
/*              and first K+1 columns of W */

#line 741 "slasyf_rook.f"
		sswap_(&k, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 742 "slasyf_rook.f"
		sswap_(&kk, &w[k + w_dim1], ldw, &w[p + w_dim1], ldw);
#line 743 "slasyf_rook.f"
	    }

/*           Updated column KP is already stored in column KK of W */

#line 747 "slasyf_rook.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

#line 751 "slasyf_rook.f"
		a[kp + k * a_dim1] = a[kk + k * a_dim1];
#line 752 "slasyf_rook.f"
		i__1 = kp - k - 1;
#line 752 "slasyf_rook.f"
		scopy_(&i__1, &a[k + 1 + kk * a_dim1], &c__1, &a[kp + (k + 1) 
			* a_dim1], lda);
#line 753 "slasyf_rook.f"
		i__1 = *n - kp + 1;
#line 753 "slasyf_rook.f"
		scopy_(&i__1, &a[kp + kk * a_dim1], &c__1, &a[kp + kp * 
			a_dim1], &c__1);

/*              Interchange rows KK and KP in first KK columns of A and W */

#line 757 "slasyf_rook.f"
		sswap_(&kk, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 758 "slasyf_rook.f"
		sswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 759 "slasyf_rook.f"
	    }

#line 761 "slasyf_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of W now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

/*              Store L(k) in column k of A */

#line 771 "slasyf_rook.f"
		i__1 = *n - k + 1;
#line 771 "slasyf_rook.f"
		scopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 772 "slasyf_rook.f"
		if (k < *n) {
#line 773 "slasyf_rook.f"
		    if ((d__1 = a[k + k * a_dim1], abs(d__1)) >= sfmin) {
#line 774 "slasyf_rook.f"
			r1 = 1. / a[k + k * a_dim1];
#line 775 "slasyf_rook.f"
			i__1 = *n - k;
#line 775 "slasyf_rook.f"
			sscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 776 "slasyf_rook.f"
		    } else if (a[k + k * a_dim1] != 0.) {
#line 777 "slasyf_rook.f"
			i__1 = *n;
#line 777 "slasyf_rook.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 778 "slasyf_rook.f"
			    a[ii + k * a_dim1] /= a[k + k * a_dim1];
#line 779 "slasyf_rook.f"
/* L74: */
#line 779 "slasyf_rook.f"
			}
#line 780 "slasyf_rook.f"
		    }
#line 781 "slasyf_rook.f"
		}

#line 783 "slasyf_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

#line 792 "slasyf_rook.f"
		if (k < *n - 1) {

/*                 Store L(k) and L(k+1) in columns k and k+1 of A */

#line 796 "slasyf_rook.f"
		    d21 = w[k + 1 + k * w_dim1];
#line 797 "slasyf_rook.f"
		    d11 = w[k + 1 + (k + 1) * w_dim1] / d21;
#line 798 "slasyf_rook.f"
		    d22 = w[k + k * w_dim1] / d21;
#line 799 "slasyf_rook.f"
		    t = 1. / (d11 * d22 - 1.);
#line 800 "slasyf_rook.f"
		    i__1 = *n;
#line 800 "slasyf_rook.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 801 "slasyf_rook.f"
			a[j + k * a_dim1] = t * ((d11 * w[j + k * w_dim1] - w[
				j + (k + 1) * w_dim1]) / d21);
#line 803 "slasyf_rook.f"
			a[j + (k + 1) * a_dim1] = t * ((d22 * w[j + (k + 1) * 
				w_dim1] - w[j + k * w_dim1]) / d21);
#line 805 "slasyf_rook.f"
/* L80: */
#line 805 "slasyf_rook.f"
		    }
#line 806 "slasyf_rook.f"
		}

/*              Copy D(k) to A */

#line 810 "slasyf_rook.f"
		a[k + k * a_dim1] = w[k + k * w_dim1];
#line 811 "slasyf_rook.f"
		a[k + 1 + k * a_dim1] = w[k + 1 + k * w_dim1];
#line 812 "slasyf_rook.f"
		a[k + 1 + (k + 1) * a_dim1] = w[k + 1 + (k + 1) * w_dim1];
#line 813 "slasyf_rook.f"
	    }
#line 814 "slasyf_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 818 "slasyf_rook.f"
	if (kstep == 1) {
#line 819 "slasyf_rook.f"
	    ipiv[k] = kp;
#line 820 "slasyf_rook.f"
	} else {
#line 821 "slasyf_rook.f"
	    ipiv[k] = -p;
#line 822 "slasyf_rook.f"
	    ipiv[k + 1] = -kp;
#line 823 "slasyf_rook.f"
	}

/*        Increase K and return to the start of the main loop */

#line 827 "slasyf_rook.f"
	k += kstep;
#line 828 "slasyf_rook.f"
	goto L70;

#line 830 "slasyf_rook.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T */

/*        computing blocks of NB columns at a time */

#line 838 "slasyf_rook.f"
	i__1 = *n;
#line 838 "slasyf_rook.f"
	i__2 = *nb;
#line 838 "slasyf_rook.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 839 "slasyf_rook.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 839 "slasyf_rook.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 843 "slasyf_rook.f"
	    i__3 = j + jb - 1;
#line 843 "slasyf_rook.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 844 "slasyf_rook.f"
		i__4 = j + jb - jj;
#line 844 "slasyf_rook.f"
		i__5 = k - 1;
#line 844 "slasyf_rook.f"
		sgemv_("No transpose", &i__4, &i__5, &c_b9, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b10, &a[jj + jj * 
			a_dim1], &c__1, (ftnlen)12);
#line 847 "slasyf_rook.f"
/* L100: */
#line 847 "slasyf_rook.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 851 "slasyf_rook.f"
	    if (j + jb <= *n) {
#line 851 "slasyf_rook.f"
		i__3 = *n - j - jb + 1;
#line 851 "slasyf_rook.f"
		i__4 = k - 1;
#line 851 "slasyf_rook.f"
		sgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &c_b9, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b10,
			 &a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 851 "slasyf_rook.f"
	    }
#line 855 "slasyf_rook.f"
/* L110: */
#line 855 "slasyf_rook.f"
	}

/*        Put L21 in standard form by partially undoing the interchanges */
/*        in columns 1:k-1 */

#line 860 "slasyf_rook.f"
	j = k - 1;
#line 861 "slasyf_rook.f"
L120:

#line 863 "slasyf_rook.f"
	kstep = 1;
#line 864 "slasyf_rook.f"
	jp1 = 1;
#line 865 "slasyf_rook.f"
	jj = j;
#line 866 "slasyf_rook.f"
	jp2 = ipiv[j];
#line 867 "slasyf_rook.f"
	if (jp2 < 0) {
#line 868 "slasyf_rook.f"
	    jp2 = -jp2;
#line 869 "slasyf_rook.f"
	    --j;
#line 870 "slasyf_rook.f"
	    jp1 = -ipiv[j];
#line 871 "slasyf_rook.f"
	    kstep = 2;
#line 872 "slasyf_rook.f"
	}

#line 874 "slasyf_rook.f"
	--j;
#line 875 "slasyf_rook.f"
	if (jp2 != jj && j >= 1) {
#line 875 "slasyf_rook.f"
	    sswap_(&j, &a[jp2 + a_dim1], lda, &a[jj + a_dim1], lda);
#line 875 "slasyf_rook.f"
	}
#line 877 "slasyf_rook.f"
	jj = j + 1;
#line 878 "slasyf_rook.f"
	if (jp1 != jj && kstep == 2) {
#line 878 "slasyf_rook.f"
	    sswap_(&j, &a[jp1 + a_dim1], lda, &a[jj + a_dim1], lda);
#line 878 "slasyf_rook.f"
	}
#line 880 "slasyf_rook.f"
	if (j >= 1) {
#line 880 "slasyf_rook.f"
	    goto L120;
#line 880 "slasyf_rook.f"
	}

/*        Set KB to the number of columns factorized */

#line 885 "slasyf_rook.f"
	*kb = k - 1;

#line 887 "slasyf_rook.f"
    }
#line 888 "slasyf_rook.f"
    return 0;

/*     End of SLASYF_ROOK */

} /* slasyf_rook__ */

