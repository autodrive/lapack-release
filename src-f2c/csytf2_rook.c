#line 1 "csytf2_rook.f"
/* csytf2_rook.f -- translated by f2c (version 20100827).
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

#line 1 "csytf2_rook.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CSYTF2_ROOK computes the factorization of a complex symmetric indefinite matrix using the bound
ed Bunch-Kaufman ("rook") diagonal pivoting method (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTF2_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytf2_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytf2_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytf2_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTF2_ROOK( UPLO, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTF2_ROOK computes the factorization of a complex symmetric matrix A */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \ingroup complexSYcomputational */

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
/* Subroutine */ int csytf2_rook__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info, ftnlen uplo_len)
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
    extern /* Subroutine */ int csyr_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal sfmin;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer itemp, kstep;
    static doublereal stemp;
    static logical upper;
    static doublereal absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 250 "csytf2_rook.f"
    /* Parameter adjustments */
#line 250 "csytf2_rook.f"
    a_dim1 = *lda;
#line 250 "csytf2_rook.f"
    a_offset = 1 + a_dim1;
#line 250 "csytf2_rook.f"
    a -= a_offset;
#line 250 "csytf2_rook.f"
    --ipiv;
#line 250 "csytf2_rook.f"

#line 250 "csytf2_rook.f"
    /* Function Body */
#line 250 "csytf2_rook.f"
    *info = 0;
#line 251 "csytf2_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 252 "csytf2_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 253 "csytf2_rook.f"
	*info = -1;
#line 254 "csytf2_rook.f"
    } else if (*n < 0) {
#line 255 "csytf2_rook.f"
	*info = -2;
#line 256 "csytf2_rook.f"
    } else if (*lda < max(1,*n)) {
#line 257 "csytf2_rook.f"
	*info = -4;
#line 258 "csytf2_rook.f"
    }
#line 259 "csytf2_rook.f"
    if (*info != 0) {
#line 260 "csytf2_rook.f"
	i__1 = -(*info);
#line 260 "csytf2_rook.f"
	xerbla_("CSYTF2_ROOK", &i__1, (ftnlen)11);
#line 261 "csytf2_rook.f"
	return 0;
#line 262 "csytf2_rook.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 266 "csytf2_rook.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 270 "csytf2_rook.f"
    sfmin = slamch_("S", (ftnlen)1);

#line 272 "csytf2_rook.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 279 "csytf2_rook.f"
	k = *n;
#line 280 "csytf2_rook.f"
L10:

/*        If K < 1, exit from loop */

#line 284 "csytf2_rook.f"
	if (k < 1) {
#line 284 "csytf2_rook.f"
	    goto L70;
#line 284 "csytf2_rook.f"
	}
#line 286 "csytf2_rook.f"
	kstep = 1;
#line 287 "csytf2_rook.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 292 "csytf2_rook.f"
	i__1 = k + k * a_dim1;
#line 292 "csytf2_rook.f"
	absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * 
		a_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 298 "csytf2_rook.f"
	if (k > 1) {
#line 299 "csytf2_rook.f"
	    i__1 = k - 1;
#line 299 "csytf2_rook.f"
	    imax = icamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 300 "csytf2_rook.f"
	    i__1 = imax + k * a_dim1;
#line 300 "csytf2_rook.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 301 "csytf2_rook.f"
	} else {
#line 302 "csytf2_rook.f"
	    colmax = 0.;
#line 303 "csytf2_rook.f"
	}

#line 305 "csytf2_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 309 "csytf2_rook.f"
	    if (*info == 0) {
#line 309 "csytf2_rook.f"
		*info = k;
#line 309 "csytf2_rook.f"
	    }
#line 311 "csytf2_rook.f"
	    kp = k;
#line 312 "csytf2_rook.f"
	} else {

/*           Test for interchange */

/*           Equivalent to testing for (used to handle NaN and Inf) */
/*           ABSAKK.GE.ALPHA*COLMAX */

#line 319 "csytf2_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, */
/*              use 1-by-1 pivot block */

#line 324 "csytf2_rook.f"
		kp = k;
#line 325 "csytf2_rook.f"
	    } else {

#line 327 "csytf2_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 331 "csytf2_rook.f"
L12:

/*                 Begin pivot search loop body */

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 339 "csytf2_rook.f"
		if (imax != k) {
#line 340 "csytf2_rook.f"
		    i__1 = k - imax;
#line 340 "csytf2_rook.f"
		    jmax = imax + icamax_(&i__1, &a[imax + (imax + 1) * 
			    a_dim1], lda);
#line 342 "csytf2_rook.f"
		    i__1 = imax + jmax * a_dim1;
#line 342 "csytf2_rook.f"
		    rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    a[imax + jmax * a_dim1]), abs(d__2));
#line 343 "csytf2_rook.f"
		} else {
#line 344 "csytf2_rook.f"
		    rowmax = 0.;
#line 345 "csytf2_rook.f"
		}

#line 347 "csytf2_rook.f"
		if (imax > 1) {
#line 348 "csytf2_rook.f"
		    i__1 = imax - 1;
#line 348 "csytf2_rook.f"
		    itemp = icamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
#line 349 "csytf2_rook.f"
		    i__1 = itemp + imax * a_dim1;
#line 349 "csytf2_rook.f"
		    stemp = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    itemp + imax * a_dim1]), abs(d__2));
#line 350 "csytf2_rook.f"
		    if (stemp > rowmax) {
#line 351 "csytf2_rook.f"
			rowmax = stemp;
#line 352 "csytf2_rook.f"
			jmax = itemp;
#line 353 "csytf2_rook.f"
		    }
#line 354 "csytf2_rook.f"
		}

/*                 Equivalent to testing for (used to handle NaN and Inf) */
/*                 CABS1( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */

#line 359 "csytf2_rook.f"
		i__1 = imax + imax * a_dim1;
#line 359 "csytf2_rook.f"
		if (! ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax 
			+ imax * a_dim1]), abs(d__2)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 365 "csytf2_rook.f"
		    kp = imax;
#line 366 "csytf2_rook.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX .EQ. COLMAX, */
/*                 used to handle NaN and Inf */

#line 371 "csytf2_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 376 "csytf2_rook.f"
		    kp = imax;
#line 377 "csytf2_rook.f"
		    kstep = 2;
#line 378 "csytf2_rook.f"
		    done = TRUE_;
#line 379 "csytf2_rook.f"
		} else {

/*                    Pivot NOT found, set variables and repeat */

#line 383 "csytf2_rook.f"
		    p = imax;
#line 384 "csytf2_rook.f"
		    colmax = rowmax;
#line 385 "csytf2_rook.f"
		    imax = jmax;
#line 386 "csytf2_rook.f"
		}

/*                 End pivot search loop body */

#line 390 "csytf2_rook.f"
		if (! done) {
#line 390 "csytf2_rook.f"
		    goto L12;
#line 390 "csytf2_rook.f"
		}

#line 392 "csytf2_rook.f"
	    }

/*           Swap TWO rows and TWO columns */

/*           First swap */

#line 398 "csytf2_rook.f"
	    if (kstep == 2 && p != k) {

/*              Interchange rows and column K and P in the leading */
/*              submatrix A(1:k,1:k) if we have a 2-by-2 pivot */

#line 403 "csytf2_rook.f"
		if (p > 1) {
#line 403 "csytf2_rook.f"
		    i__1 = p - 1;
#line 403 "csytf2_rook.f"
		    cswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 
			    1], &c__1);
#line 403 "csytf2_rook.f"
		}
#line 405 "csytf2_rook.f"
		if (p < k - 1) {
#line 405 "csytf2_rook.f"
		    i__1 = k - p - 1;
#line 405 "csytf2_rook.f"
		    cswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 
			    1) * a_dim1], lda);
#line 405 "csytf2_rook.f"
		}
#line 408 "csytf2_rook.f"
		i__1 = k + k * a_dim1;
#line 408 "csytf2_rook.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 409 "csytf2_rook.f"
		i__1 = k + k * a_dim1;
#line 409 "csytf2_rook.f"
		i__2 = p + p * a_dim1;
#line 409 "csytf2_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 410 "csytf2_rook.f"
		i__1 = p + p * a_dim1;
#line 410 "csytf2_rook.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 411 "csytf2_rook.f"
	    }

/*           Second swap */

#line 415 "csytf2_rook.f"
	    kk = k - kstep + 1;
#line 416 "csytf2_rook.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 421 "csytf2_rook.f"
		if (kp > 1) {
#line 421 "csytf2_rook.f"
		    i__1 = kp - 1;
#line 421 "csytf2_rook.f"
		    cswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 421 "csytf2_rook.f"
		}
#line 423 "csytf2_rook.f"
		if (kk > 1 && kp < kk - 1) {
#line 423 "csytf2_rook.f"
		    i__1 = kk - kp - 1;
#line 423 "csytf2_rook.f"
		    cswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (
			    kp + 1) * a_dim1], lda);
#line 423 "csytf2_rook.f"
		}
#line 426 "csytf2_rook.f"
		i__1 = kk + kk * a_dim1;
#line 426 "csytf2_rook.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 427 "csytf2_rook.f"
		i__1 = kk + kk * a_dim1;
#line 427 "csytf2_rook.f"
		i__2 = kp + kp * a_dim1;
#line 427 "csytf2_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 428 "csytf2_rook.f"
		i__1 = kp + kp * a_dim1;
#line 428 "csytf2_rook.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 429 "csytf2_rook.f"
		if (kstep == 2) {
#line 430 "csytf2_rook.f"
		    i__1 = k - 1 + k * a_dim1;
#line 430 "csytf2_rook.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 431 "csytf2_rook.f"
		    i__1 = k - 1 + k * a_dim1;
#line 431 "csytf2_rook.f"
		    i__2 = kp + k * a_dim1;
#line 431 "csytf2_rook.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 432 "csytf2_rook.f"
		    i__1 = kp + k * a_dim1;
#line 432 "csytf2_rook.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 433 "csytf2_rook.f"
		}
#line 434 "csytf2_rook.f"
	    }

/*           Update the leading submatrix */

#line 438 "csytf2_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

#line 446 "csytf2_rook.f"
		if (k > 1) {

/*                 Perform a rank-1 update of A(1:k-1,1:k-1) and */
/*                 store U(k) in column k */

#line 451 "csytf2_rook.f"
		    i__1 = k + k * a_dim1;
#line 451 "csytf2_rook.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + 
			    k * a_dim1]), abs(d__2)) >= sfmin) {

/*                    Perform a rank-1 update of A(1:k-1,1:k-1) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*1/D(k)*W(k)**T */

#line 457 "csytf2_rook.f"
			z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 457 "csytf2_rook.f"
			d11.r = z__1.r, d11.i = z__1.i;
#line 458 "csytf2_rook.f"
			i__1 = k - 1;
#line 458 "csytf2_rook.f"
			z__1.r = -d11.r, z__1.i = -d11.i;
#line 458 "csytf2_rook.f"
			csyr_(uplo, &i__1, &z__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);

/*                    Store U(k) in column k */

#line 462 "csytf2_rook.f"
			i__1 = k - 1;
#line 462 "csytf2_rook.f"
			cscal_(&i__1, &d11, &a[k * a_dim1 + 1], &c__1);
#line 463 "csytf2_rook.f"
		    } else {

/*                    Store L(k) in column K */

#line 467 "csytf2_rook.f"
			i__1 = k + k * a_dim1;
#line 467 "csytf2_rook.f"
			d11.r = a[i__1].r, d11.i = a[i__1].i;
#line 468 "csytf2_rook.f"
			i__1 = k - 1;
#line 468 "csytf2_rook.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 469 "csytf2_rook.f"
			    i__2 = ii + k * a_dim1;
#line 469 "csytf2_rook.f"
			    z_div(&z__1, &a[ii + k * a_dim1], &d11);
#line 469 "csytf2_rook.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 470 "csytf2_rook.f"
/* L16: */
#line 470 "csytf2_rook.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 477 "csytf2_rook.f"
			i__1 = k - 1;
#line 477 "csytf2_rook.f"
			z__1.r = -d11.r, z__1.i = -d11.i;
#line 477 "csytf2_rook.f"
			csyr_(uplo, &i__1, &z__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);
#line 478 "csytf2_rook.f"
		    }
#line 479 "csytf2_rook.f"
		}

#line 481 "csytf2_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 497 "csytf2_rook.f"
		if (k > 2) {

#line 499 "csytf2_rook.f"
		    i__1 = k - 1 + k * a_dim1;
#line 499 "csytf2_rook.f"
		    d12.r = a[i__1].r, d12.i = a[i__1].i;
#line 500 "csytf2_rook.f"
		    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &d12);
#line 500 "csytf2_rook.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 501 "csytf2_rook.f"
		    z_div(&z__1, &a[k + k * a_dim1], &d12);
#line 501 "csytf2_rook.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 502 "csytf2_rook.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 502 "csytf2_rook.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 502 "csytf2_rook.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 502 "csytf2_rook.f"
		    t.r = z__1.r, t.i = z__1.i;

#line 504 "csytf2_rook.f"
		    for (j = k - 2; j >= 1; --j) {

#line 506 "csytf2_rook.f"
			i__1 = j + (k - 1) * a_dim1;
#line 506 "csytf2_rook.f"
			z__3.r = d11.r * a[i__1].r - d11.i * a[i__1].i, 
				z__3.i = d11.r * a[i__1].i + d11.i * a[i__1]
				.r;
#line 506 "csytf2_rook.f"
			i__2 = j + k * a_dim1;
#line 506 "csytf2_rook.f"
			z__2.r = z__3.r - a[i__2].r, z__2.i = z__3.i - a[i__2]
				.i;
#line 506 "csytf2_rook.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 506 "csytf2_rook.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 507 "csytf2_rook.f"
			i__1 = j + k * a_dim1;
#line 507 "csytf2_rook.f"
			z__3.r = d22.r * a[i__1].r - d22.i * a[i__1].i, 
				z__3.i = d22.r * a[i__1].i + d22.i * a[i__1]
				.r;
#line 507 "csytf2_rook.f"
			i__2 = j + (k - 1) * a_dim1;
#line 507 "csytf2_rook.f"
			z__2.r = z__3.r - a[i__2].r, z__2.i = z__3.i - a[i__2]
				.i;
#line 507 "csytf2_rook.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 507 "csytf2_rook.f"
			wk.r = z__1.r, wk.i = z__1.i;

#line 509 "csytf2_rook.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 510 "csytf2_rook.f"
			    i__1 = i__ + j * a_dim1;
#line 510 "csytf2_rook.f"
			    i__2 = i__ + j * a_dim1;
#line 510 "csytf2_rook.f"
			    z_div(&z__4, &a[i__ + k * a_dim1], &d12);
#line 510 "csytf2_rook.f"
			    z__3.r = z__4.r * wk.r - z__4.i * wk.i, z__3.i = 
				    z__4.r * wk.i + z__4.i * wk.r;
#line 510 "csytf2_rook.f"
			    z__2.r = a[i__2].r - z__3.r, z__2.i = a[i__2].i - 
				    z__3.i;
#line 510 "csytf2_rook.f"
			    z_div(&z__6, &a[i__ + (k - 1) * a_dim1], &d12);
#line 510 "csytf2_rook.f"
			    z__5.r = z__6.r * wkm1.r - z__6.i * wkm1.i, 
				    z__5.i = z__6.r * wkm1.i + z__6.i * 
				    wkm1.r;
#line 510 "csytf2_rook.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 510 "csytf2_rook.f"
			    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 512 "csytf2_rook.f"
/* L20: */
#line 512 "csytf2_rook.f"
			}

/*                    Store U(k) and U(k-1) in cols k and k-1 for row J */

#line 516 "csytf2_rook.f"
			i__1 = j + k * a_dim1;
#line 516 "csytf2_rook.f"
			z_div(&z__1, &wk, &d12);
#line 516 "csytf2_rook.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 517 "csytf2_rook.f"
			i__1 = j + (k - 1) * a_dim1;
#line 517 "csytf2_rook.f"
			z_div(&z__1, &wkm1, &d12);
#line 517 "csytf2_rook.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 519 "csytf2_rook.f"
/* L30: */
#line 519 "csytf2_rook.f"
		    }

#line 521 "csytf2_rook.f"
		}

#line 523 "csytf2_rook.f"
	    }
#line 524 "csytf2_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 528 "csytf2_rook.f"
	if (kstep == 1) {
#line 529 "csytf2_rook.f"
	    ipiv[k] = kp;
#line 530 "csytf2_rook.f"
	} else {
#line 531 "csytf2_rook.f"
	    ipiv[k] = -p;
#line 532 "csytf2_rook.f"
	    ipiv[k - 1] = -kp;
#line 533 "csytf2_rook.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 537 "csytf2_rook.f"
	k -= kstep;
#line 538 "csytf2_rook.f"
	goto L10;

#line 540 "csytf2_rook.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 547 "csytf2_rook.f"
	k = 1;
#line 548 "csytf2_rook.f"
L40:

/*        If K > N, exit from loop */

#line 552 "csytf2_rook.f"
	if (k > *n) {
#line 552 "csytf2_rook.f"
	    goto L70;
#line 552 "csytf2_rook.f"
	}
#line 554 "csytf2_rook.f"
	kstep = 1;
#line 555 "csytf2_rook.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 560 "csytf2_rook.f"
	i__1 = k + k * a_dim1;
#line 560 "csytf2_rook.f"
	absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * 
		a_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 566 "csytf2_rook.f"
	if (k < *n) {
#line 567 "csytf2_rook.f"
	    i__1 = *n - k;
#line 567 "csytf2_rook.f"
	    imax = k + icamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 568 "csytf2_rook.f"
	    i__1 = imax + k * a_dim1;
#line 568 "csytf2_rook.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 569 "csytf2_rook.f"
	} else {
#line 570 "csytf2_rook.f"
	    colmax = 0.;
#line 571 "csytf2_rook.f"
	}

#line 573 "csytf2_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 577 "csytf2_rook.f"
	    if (*info == 0) {
#line 577 "csytf2_rook.f"
		*info = k;
#line 577 "csytf2_rook.f"
	    }
#line 579 "csytf2_rook.f"
	    kp = k;
#line 580 "csytf2_rook.f"
	} else {

/*           Test for interchange */

/*           Equivalent to testing for (used to handle NaN and Inf) */
/*           ABSAKK.GE.ALPHA*COLMAX */

#line 587 "csytf2_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 591 "csytf2_rook.f"
		kp = k;
#line 592 "csytf2_rook.f"
	    } else {

#line 594 "csytf2_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 598 "csytf2_rook.f"
L42:

/*                 Begin pivot search loop body */

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 606 "csytf2_rook.f"
		if (imax != k) {
#line 607 "csytf2_rook.f"
		    i__1 = imax - k;
#line 607 "csytf2_rook.f"
		    jmax = k - 1 + icamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 608 "csytf2_rook.f"
		    i__1 = imax + jmax * a_dim1;
#line 608 "csytf2_rook.f"
		    rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    a[imax + jmax * a_dim1]), abs(d__2));
#line 609 "csytf2_rook.f"
		} else {
#line 610 "csytf2_rook.f"
		    rowmax = 0.;
#line 611 "csytf2_rook.f"
		}

#line 613 "csytf2_rook.f"
		if (imax < *n) {
#line 614 "csytf2_rook.f"
		    i__1 = *n - imax;
#line 614 "csytf2_rook.f"
		    itemp = imax + icamax_(&i__1, &a[imax + 1 + imax * a_dim1]
			    , &c__1);
#line 616 "csytf2_rook.f"
		    i__1 = itemp + imax * a_dim1;
#line 616 "csytf2_rook.f"
		    stemp = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    itemp + imax * a_dim1]), abs(d__2));
#line 617 "csytf2_rook.f"
		    if (stemp > rowmax) {
#line 618 "csytf2_rook.f"
			rowmax = stemp;
#line 619 "csytf2_rook.f"
			jmax = itemp;
#line 620 "csytf2_rook.f"
		    }
#line 621 "csytf2_rook.f"
		}

/*                 Equivalent to testing for (used to handle NaN and Inf) */
/*                 CABS1( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */

#line 626 "csytf2_rook.f"
		i__1 = imax + imax * a_dim1;
#line 626 "csytf2_rook.f"
		if (! ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax 
			+ imax * a_dim1]), abs(d__2)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 632 "csytf2_rook.f"
		    kp = imax;
#line 633 "csytf2_rook.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX .EQ. COLMAX, */
/*                 used to handle NaN and Inf */

#line 638 "csytf2_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 643 "csytf2_rook.f"
		    kp = imax;
#line 644 "csytf2_rook.f"
		    kstep = 2;
#line 645 "csytf2_rook.f"
		    done = TRUE_;
#line 646 "csytf2_rook.f"
		} else {

/*                    Pivot NOT found, set variables and repeat */

#line 650 "csytf2_rook.f"
		    p = imax;
#line 651 "csytf2_rook.f"
		    colmax = rowmax;
#line 652 "csytf2_rook.f"
		    imax = jmax;
#line 653 "csytf2_rook.f"
		}

/*                 End pivot search loop body */

#line 657 "csytf2_rook.f"
		if (! done) {
#line 657 "csytf2_rook.f"
		    goto L42;
#line 657 "csytf2_rook.f"
		}

#line 659 "csytf2_rook.f"
	    }

/*           Swap TWO rows and TWO columns */

/*           First swap */

#line 665 "csytf2_rook.f"
	    if (kstep == 2 && p != k) {

/*              Interchange rows and column K and P in the trailing */
/*              submatrix A(k:n,k:n) if we have a 2-by-2 pivot */

#line 670 "csytf2_rook.f"
		if (p < *n) {
#line 670 "csytf2_rook.f"
		    i__1 = *n - p;
#line 670 "csytf2_rook.f"
		    cswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p 
			    * a_dim1], &c__1);
#line 670 "csytf2_rook.f"
		}
#line 672 "csytf2_rook.f"
		if (p > k + 1) {
#line 672 "csytf2_rook.f"
		    i__1 = p - k - 1;
#line 672 "csytf2_rook.f"
		    cswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[p + (k + 
			    1) * a_dim1], lda);
#line 672 "csytf2_rook.f"
		}
#line 674 "csytf2_rook.f"
		i__1 = k + k * a_dim1;
#line 674 "csytf2_rook.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 675 "csytf2_rook.f"
		i__1 = k + k * a_dim1;
#line 675 "csytf2_rook.f"
		i__2 = p + p * a_dim1;
#line 675 "csytf2_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 676 "csytf2_rook.f"
		i__1 = p + p * a_dim1;
#line 676 "csytf2_rook.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 677 "csytf2_rook.f"
	    }

/*           Second swap */

#line 681 "csytf2_rook.f"
	    kk = k + kstep - 1;
#line 682 "csytf2_rook.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 687 "csytf2_rook.f"
		if (kp < *n) {
#line 687 "csytf2_rook.f"
		    i__1 = *n - kp;
#line 687 "csytf2_rook.f"
		    cswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 687 "csytf2_rook.f"
		}
#line 689 "csytf2_rook.f"
		if (kk < *n && kp > kk + 1) {
#line 689 "csytf2_rook.f"
		    i__1 = kp - kk - 1;
#line 689 "csytf2_rook.f"
		    cswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (
			    kk + 1) * a_dim1], lda);
#line 689 "csytf2_rook.f"
		}
#line 692 "csytf2_rook.f"
		i__1 = kk + kk * a_dim1;
#line 692 "csytf2_rook.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 693 "csytf2_rook.f"
		i__1 = kk + kk * a_dim1;
#line 693 "csytf2_rook.f"
		i__2 = kp + kp * a_dim1;
#line 693 "csytf2_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 694 "csytf2_rook.f"
		i__1 = kp + kp * a_dim1;
#line 694 "csytf2_rook.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 695 "csytf2_rook.f"
		if (kstep == 2) {
#line 696 "csytf2_rook.f"
		    i__1 = k + 1 + k * a_dim1;
#line 696 "csytf2_rook.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 697 "csytf2_rook.f"
		    i__1 = k + 1 + k * a_dim1;
#line 697 "csytf2_rook.f"
		    i__2 = kp + k * a_dim1;
#line 697 "csytf2_rook.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 698 "csytf2_rook.f"
		    i__1 = kp + k * a_dim1;
#line 698 "csytf2_rook.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 699 "csytf2_rook.f"
		}
#line 700 "csytf2_rook.f"
	    }

/*           Update the trailing submatrix */

#line 704 "csytf2_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 712 "csytf2_rook.f"
		if (k < *n) {

/*              Perform a rank-1 update of A(k+1:n,k+1:n) and */
/*              store L(k) in column k */

#line 717 "csytf2_rook.f"
		    i__1 = k + k * a_dim1;
#line 717 "csytf2_rook.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + 
			    k * a_dim1]), abs(d__2)) >= sfmin) {

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */

#line 723 "csytf2_rook.f"
			z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 723 "csytf2_rook.f"
			d11.r = z__1.r, d11.i = z__1.i;
#line 724 "csytf2_rook.f"
			i__1 = *n - k;
#line 724 "csytf2_rook.f"
			z__1.r = -d11.r, z__1.i = -d11.i;
#line 724 "csytf2_rook.f"
			csyr_(uplo, &i__1, &z__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);

/*                    Store L(k) in column k */

#line 729 "csytf2_rook.f"
			i__1 = *n - k;
#line 729 "csytf2_rook.f"
			cscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
#line 730 "csytf2_rook.f"
		    } else {

/*                    Store L(k) in column k */

#line 734 "csytf2_rook.f"
			i__1 = k + k * a_dim1;
#line 734 "csytf2_rook.f"
			d11.r = a[i__1].r, d11.i = a[i__1].i;
#line 735 "csytf2_rook.f"
			i__1 = *n;
#line 735 "csytf2_rook.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 736 "csytf2_rook.f"
			    i__2 = ii + k * a_dim1;
#line 736 "csytf2_rook.f"
			    z_div(&z__1, &a[ii + k * a_dim1], &d11);
#line 736 "csytf2_rook.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 737 "csytf2_rook.f"
/* L46: */
#line 737 "csytf2_rook.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 744 "csytf2_rook.f"
			i__1 = *n - k;
#line 744 "csytf2_rook.f"
			z__1.r = -d11.r, z__1.i = -d11.i;
#line 744 "csytf2_rook.f"
			csyr_(uplo, &i__1, &z__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);
#line 746 "csytf2_rook.f"
		    }
#line 747 "csytf2_rook.f"
		}

#line 749 "csytf2_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */


/*              Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T */
/*                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 766 "csytf2_rook.f"
		if (k < *n - 1) {

#line 768 "csytf2_rook.f"
		    i__1 = k + 1 + k * a_dim1;
#line 768 "csytf2_rook.f"
		    d21.r = a[i__1].r, d21.i = a[i__1].i;
#line 769 "csytf2_rook.f"
		    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &d21);
#line 769 "csytf2_rook.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 770 "csytf2_rook.f"
		    z_div(&z__1, &a[k + k * a_dim1], &d21);
#line 770 "csytf2_rook.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 771 "csytf2_rook.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 771 "csytf2_rook.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 771 "csytf2_rook.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 771 "csytf2_rook.f"
		    t.r = z__1.r, t.i = z__1.i;

#line 773 "csytf2_rook.f"
		    i__1 = *n;
#line 773 "csytf2_rook.f"
		    for (j = k + 2; j <= i__1; ++j) {

/*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */

#line 777 "csytf2_rook.f"
			i__2 = j + k * a_dim1;
#line 777 "csytf2_rook.f"
			z__3.r = d11.r * a[i__2].r - d11.i * a[i__2].i, 
				z__3.i = d11.r * a[i__2].i + d11.i * a[i__2]
				.r;
#line 777 "csytf2_rook.f"
			i__3 = j + (k + 1) * a_dim1;
#line 777 "csytf2_rook.f"
			z__2.r = z__3.r - a[i__3].r, z__2.i = z__3.i - a[i__3]
				.i;
#line 777 "csytf2_rook.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 777 "csytf2_rook.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 778 "csytf2_rook.f"
			i__2 = j + (k + 1) * a_dim1;
#line 778 "csytf2_rook.f"
			z__3.r = d22.r * a[i__2].r - d22.i * a[i__2].i, 
				z__3.i = d22.r * a[i__2].i + d22.i * a[i__2]
				.r;
#line 778 "csytf2_rook.f"
			i__3 = j + k * a_dim1;
#line 778 "csytf2_rook.f"
			z__2.r = z__3.r - a[i__3].r, z__2.i = z__3.i - a[i__3]
				.i;
#line 778 "csytf2_rook.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 778 "csytf2_rook.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;

/*                    Perform a rank-2 update of A(k+2:n,k+2:n) */

#line 782 "csytf2_rook.f"
			i__2 = *n;
#line 782 "csytf2_rook.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 783 "csytf2_rook.f"
			    i__3 = i__ + j * a_dim1;
#line 783 "csytf2_rook.f"
			    i__4 = i__ + j * a_dim1;
#line 783 "csytf2_rook.f"
			    z_div(&z__4, &a[i__ + k * a_dim1], &d21);
#line 783 "csytf2_rook.f"
			    z__3.r = z__4.r * wk.r - z__4.i * wk.i, z__3.i = 
				    z__4.r * wk.i + z__4.i * wk.r;
#line 783 "csytf2_rook.f"
			    z__2.r = a[i__4].r - z__3.r, z__2.i = a[i__4].i - 
				    z__3.i;
#line 783 "csytf2_rook.f"
			    z_div(&z__6, &a[i__ + (k + 1) * a_dim1], &d21);
#line 783 "csytf2_rook.f"
			    z__5.r = z__6.r * wkp1.r - z__6.i * wkp1.i, 
				    z__5.i = z__6.r * wkp1.i + z__6.i * 
				    wkp1.r;
#line 783 "csytf2_rook.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 783 "csytf2_rook.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 785 "csytf2_rook.f"
/* L50: */
#line 785 "csytf2_rook.f"
			}

/*                    Store L(k) and L(k+1) in cols k and k+1 for row J */

#line 789 "csytf2_rook.f"
			i__2 = j + k * a_dim1;
#line 789 "csytf2_rook.f"
			z_div(&z__1, &wk, &d21);
#line 789 "csytf2_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 790 "csytf2_rook.f"
			i__2 = j + (k + 1) * a_dim1;
#line 790 "csytf2_rook.f"
			z_div(&z__1, &wkp1, &d21);
#line 790 "csytf2_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;

#line 792 "csytf2_rook.f"
/* L60: */
#line 792 "csytf2_rook.f"
		    }

#line 794 "csytf2_rook.f"
		}

#line 796 "csytf2_rook.f"
	    }
#line 797 "csytf2_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 801 "csytf2_rook.f"
	if (kstep == 1) {
#line 802 "csytf2_rook.f"
	    ipiv[k] = kp;
#line 803 "csytf2_rook.f"
	} else {
#line 804 "csytf2_rook.f"
	    ipiv[k] = -p;
#line 805 "csytf2_rook.f"
	    ipiv[k + 1] = -kp;
#line 806 "csytf2_rook.f"
	}

/*        Increase K and return to the start of the main loop */

#line 810 "csytf2_rook.f"
	k += kstep;
#line 811 "csytf2_rook.f"
	goto L40;

#line 813 "csytf2_rook.f"
    }

#line 815 "csytf2_rook.f"
L70:

#line 817 "csytf2_rook.f"
    return 0;

/*     End of CSYTF2_ROOK */

} /* csytf2_rook__ */

