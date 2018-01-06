#line 1 "chetf2_rook.f"
/* chetf2_rook.f -- translated by f2c (version 20100827).
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

#line 1 "chetf2_rook.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHETF2_ROOK computes the factorization of a complex Hermitian indefinite matrix using the bound
ed Bunch-Kaufman ("rook") diagonal pivoting method (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETF2_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetf2_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetf2_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetf2_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETF2_ROOK( UPLO, N, A, LDA, IPIV, INFO ) */

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
/* > CHETF2_ROOK computes the factorization of a complex Hermitian matrix A */
/* > using the bounded Bunch-Kaufman ("rook") diagonal pivoting method: */
/* > */
/* >    A = U*D*U**H  or  A = L*D*L**H */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, U**H is the conjugate transpose of U, and D is */
/* > Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
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
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >             interchanged and D(k,k) is a 1-by-1 diagonal block. */
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

/* > \ingroup complexHEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  If UPLO = 'U', then A = U*D*U**H, where */
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
/* >  If UPLO = 'L', then A = L*D*L**H, where */
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
/* >  November 2013,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* >  01-01-96 - Based on modifications by */
/* >    J. Lewis, Boeing Computer Services Company */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int chetf2_rook__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info, ftnlen uplo_len)
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
    extern /* Subroutine */ int cher_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static logical done;
    static integer imax, jmax;
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal sfmin;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer itemp, kstep;
    static doublereal stemp;
    static logical upper;
    extern doublereal slapy2_(doublereal *, doublereal *);
    static doublereal absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    static doublereal colmax, rowmax;


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

#line 250 "chetf2_rook.f"
    /* Parameter adjustments */
#line 250 "chetf2_rook.f"
    a_dim1 = *lda;
#line 250 "chetf2_rook.f"
    a_offset = 1 + a_dim1;
#line 250 "chetf2_rook.f"
    a -= a_offset;
#line 250 "chetf2_rook.f"
    --ipiv;
#line 250 "chetf2_rook.f"

#line 250 "chetf2_rook.f"
    /* Function Body */
#line 250 "chetf2_rook.f"
    *info = 0;
#line 251 "chetf2_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 252 "chetf2_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 253 "chetf2_rook.f"
	*info = -1;
#line 254 "chetf2_rook.f"
    } else if (*n < 0) {
#line 255 "chetf2_rook.f"
	*info = -2;
#line 256 "chetf2_rook.f"
    } else if (*lda < max(1,*n)) {
#line 257 "chetf2_rook.f"
	*info = -4;
#line 258 "chetf2_rook.f"
    }
#line 259 "chetf2_rook.f"
    if (*info != 0) {
#line 260 "chetf2_rook.f"
	i__1 = -(*info);
#line 260 "chetf2_rook.f"
	xerbla_("CHETF2_ROOK", &i__1, (ftnlen)11);
#line 261 "chetf2_rook.f"
	return 0;
#line 262 "chetf2_rook.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 266 "chetf2_rook.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 270 "chetf2_rook.f"
    sfmin = slamch_("S", (ftnlen)1);

#line 272 "chetf2_rook.f"
    if (upper) {

/*        Factorize A as U*D*U**H using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 279 "chetf2_rook.f"
	k = *n;
#line 280 "chetf2_rook.f"
L10:

/*        If K < 1, exit from loop */

#line 284 "chetf2_rook.f"
	if (k < 1) {
#line 284 "chetf2_rook.f"
	    goto L70;
#line 284 "chetf2_rook.f"
	}
#line 286 "chetf2_rook.f"
	kstep = 1;
#line 287 "chetf2_rook.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 292 "chetf2_rook.f"
	i__1 = k + k * a_dim1;
#line 292 "chetf2_rook.f"
	absakk = (d__1 = a[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 298 "chetf2_rook.f"
	if (k > 1) {
#line 299 "chetf2_rook.f"
	    i__1 = k - 1;
#line 299 "chetf2_rook.f"
	    imax = icamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 300 "chetf2_rook.f"
	    i__1 = imax + k * a_dim1;
#line 300 "chetf2_rook.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 301 "chetf2_rook.f"
	} else {
#line 302 "chetf2_rook.f"
	    colmax = 0.;
#line 303 "chetf2_rook.f"
	}

#line 305 "chetf2_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 309 "chetf2_rook.f"
	    if (*info == 0) {
#line 309 "chetf2_rook.f"
		*info = k;
#line 309 "chetf2_rook.f"
	    }
#line 311 "chetf2_rook.f"
	    kp = k;
#line 312 "chetf2_rook.f"
	    i__1 = k + k * a_dim1;
#line 312 "chetf2_rook.f"
	    i__2 = k + k * a_dim1;
#line 312 "chetf2_rook.f"
	    d__1 = a[i__2].r;
#line 312 "chetf2_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 313 "chetf2_rook.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 323 "chetf2_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 327 "chetf2_rook.f"
		kp = k;

#line 329 "chetf2_rook.f"
	    } else {

#line 331 "chetf2_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 335 "chetf2_rook.f"
L12:

/*                 BEGIN pivot search loop body */


/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 344 "chetf2_rook.f"
		if (imax != k) {
#line 345 "chetf2_rook.f"
		    i__1 = k - imax;
#line 345 "chetf2_rook.f"
		    jmax = imax + icamax_(&i__1, &a[imax + (imax + 1) * 
			    a_dim1], lda);
#line 347 "chetf2_rook.f"
		    i__1 = imax + jmax * a_dim1;
#line 347 "chetf2_rook.f"
		    rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    a[imax + jmax * a_dim1]), abs(d__2));
#line 348 "chetf2_rook.f"
		} else {
#line 349 "chetf2_rook.f"
		    rowmax = 0.;
#line 350 "chetf2_rook.f"
		}

#line 352 "chetf2_rook.f"
		if (imax > 1) {
#line 353 "chetf2_rook.f"
		    i__1 = imax - 1;
#line 353 "chetf2_rook.f"
		    itemp = icamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
#line 354 "chetf2_rook.f"
		    i__1 = itemp + imax * a_dim1;
#line 354 "chetf2_rook.f"
		    stemp = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    itemp + imax * a_dim1]), abs(d__2));
#line 355 "chetf2_rook.f"
		    if (stemp > rowmax) {
#line 356 "chetf2_rook.f"
			rowmax = stemp;
#line 357 "chetf2_rook.f"
			jmax = itemp;
#line 358 "chetf2_rook.f"
		    }
#line 359 "chetf2_rook.f"
		}

/*                 Case(2) */
/*                 Equivalent to testing for */
/*                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 366 "chetf2_rook.f"
		i__1 = imax + imax * a_dim1;
#line 366 "chetf2_rook.f"
		if (! ((d__1 = a[i__1].r, abs(d__1)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 372 "chetf2_rook.f"
		    kp = imax;
#line 373 "chetf2_rook.f"
		    done = TRUE_;

/*                 Case(3) */
/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 379 "chetf2_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 385 "chetf2_rook.f"
		    kp = imax;
#line 386 "chetf2_rook.f"
		    kstep = 2;
#line 387 "chetf2_rook.f"
		    done = TRUE_;

/*                 Case(4) */
#line 390 "chetf2_rook.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 394 "chetf2_rook.f"
		    p = imax;
#line 395 "chetf2_rook.f"
		    colmax = rowmax;
#line 396 "chetf2_rook.f"
		    imax = jmax;
#line 397 "chetf2_rook.f"
		}

/*                 END pivot search loop body */

#line 401 "chetf2_rook.f"
		if (! done) {
#line 401 "chetf2_rook.f"
		    goto L12;
#line 401 "chetf2_rook.f"
		}

#line 403 "chetf2_rook.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 411 "chetf2_rook.f"
	    kk = k - kstep + 1;

/*           For only a 2x2 pivot, interchange rows and columns K and P */
/*           in the leading submatrix A(1:k,1:k) */

#line 416 "chetf2_rook.f"
	    if (kstep == 2 && p != k) {
/*              (1) Swap columnar parts */
#line 418 "chetf2_rook.f"
		if (p > 1) {
#line 418 "chetf2_rook.f"
		    i__1 = p - 1;
#line 418 "chetf2_rook.f"
		    cswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 
			    1], &c__1);
#line 418 "chetf2_rook.f"
		}
/*              (2) Swap and conjugate middle parts */
#line 421 "chetf2_rook.f"
		i__1 = k - 1;
#line 421 "chetf2_rook.f"
		for (j = p + 1; j <= i__1; ++j) {
#line 422 "chetf2_rook.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 422 "chetf2_rook.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 423 "chetf2_rook.f"
		    i__2 = j + k * a_dim1;
#line 423 "chetf2_rook.f"
		    d_cnjg(&z__1, &a[p + j * a_dim1]);
#line 423 "chetf2_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 424 "chetf2_rook.f"
		    i__2 = p + j * a_dim1;
#line 424 "chetf2_rook.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 425 "chetf2_rook.f"
/* L14: */
#line 425 "chetf2_rook.f"
		}
/*              (3) Swap and conjugate corner elements at row-col interserction */
#line 427 "chetf2_rook.f"
		i__1 = p + k * a_dim1;
#line 427 "chetf2_rook.f"
		d_cnjg(&z__1, &a[p + k * a_dim1]);
#line 427 "chetf2_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*              (4) Swap diagonal elements at row-col intersection */
#line 429 "chetf2_rook.f"
		i__1 = k + k * a_dim1;
#line 429 "chetf2_rook.f"
		r1 = a[i__1].r;
#line 430 "chetf2_rook.f"
		i__1 = k + k * a_dim1;
#line 430 "chetf2_rook.f"
		i__2 = p + p * a_dim1;
#line 430 "chetf2_rook.f"
		d__1 = a[i__2].r;
#line 430 "chetf2_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 431 "chetf2_rook.f"
		i__1 = p + p * a_dim1;
#line 431 "chetf2_rook.f"
		a[i__1].r = r1, a[i__1].i = 0.;
#line 432 "chetf2_rook.f"
	    }

/*           For both 1x1 and 2x2 pivots, interchange rows and */
/*           columns KK and KP in the leading submatrix A(1:k,1:k) */

#line 437 "chetf2_rook.f"
	    if (kp != kk) {
/*              (1) Swap columnar parts */
#line 439 "chetf2_rook.f"
		if (kp > 1) {
#line 439 "chetf2_rook.f"
		    i__1 = kp - 1;
#line 439 "chetf2_rook.f"
		    cswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 439 "chetf2_rook.f"
		}
/*              (2) Swap and conjugate middle parts */
#line 442 "chetf2_rook.f"
		i__1 = kk - 1;
#line 442 "chetf2_rook.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 443 "chetf2_rook.f"
		    d_cnjg(&z__1, &a[j + kk * a_dim1]);
#line 443 "chetf2_rook.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 444 "chetf2_rook.f"
		    i__2 = j + kk * a_dim1;
#line 444 "chetf2_rook.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 444 "chetf2_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 445 "chetf2_rook.f"
		    i__2 = kp + j * a_dim1;
#line 445 "chetf2_rook.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 446 "chetf2_rook.f"
/* L15: */
#line 446 "chetf2_rook.f"
		}
/*              (3) Swap and conjugate corner elements at row-col interserction */
#line 448 "chetf2_rook.f"
		i__1 = kp + kk * a_dim1;
#line 448 "chetf2_rook.f"
		d_cnjg(&z__1, &a[kp + kk * a_dim1]);
#line 448 "chetf2_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*              (4) Swap diagonal elements at row-col intersection */
#line 450 "chetf2_rook.f"
		i__1 = kk + kk * a_dim1;
#line 450 "chetf2_rook.f"
		r1 = a[i__1].r;
#line 451 "chetf2_rook.f"
		i__1 = kk + kk * a_dim1;
#line 451 "chetf2_rook.f"
		i__2 = kp + kp * a_dim1;
#line 451 "chetf2_rook.f"
		d__1 = a[i__2].r;
#line 451 "chetf2_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 452 "chetf2_rook.f"
		i__1 = kp + kp * a_dim1;
#line 452 "chetf2_rook.f"
		a[i__1].r = r1, a[i__1].i = 0.;

#line 454 "chetf2_rook.f"
		if (kstep == 2) {
/*                 (*) Make sure that diagonal element of pivot is real */
#line 456 "chetf2_rook.f"
		    i__1 = k + k * a_dim1;
#line 456 "chetf2_rook.f"
		    i__2 = k + k * a_dim1;
#line 456 "chetf2_rook.f"
		    d__1 = a[i__2].r;
#line 456 "chetf2_rook.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
/*                 (5) Swap row elements */
#line 458 "chetf2_rook.f"
		    i__1 = k - 1 + k * a_dim1;
#line 458 "chetf2_rook.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 459 "chetf2_rook.f"
		    i__1 = k - 1 + k * a_dim1;
#line 459 "chetf2_rook.f"
		    i__2 = kp + k * a_dim1;
#line 459 "chetf2_rook.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 460 "chetf2_rook.f"
		    i__1 = kp + k * a_dim1;
#line 460 "chetf2_rook.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 461 "chetf2_rook.f"
		}
#line 462 "chetf2_rook.f"
	    } else {
/*              (*) Make sure that diagonal element of pivot is real */
#line 464 "chetf2_rook.f"
		i__1 = k + k * a_dim1;
#line 464 "chetf2_rook.f"
		i__2 = k + k * a_dim1;
#line 464 "chetf2_rook.f"
		d__1 = a[i__2].r;
#line 464 "chetf2_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 465 "chetf2_rook.f"
		if (kstep == 2) {
#line 465 "chetf2_rook.f"
		    i__1 = k - 1 + (k - 1) * a_dim1;
#line 465 "chetf2_rook.f"
		    i__2 = k - 1 + (k - 1) * a_dim1;
#line 465 "chetf2_rook.f"
		    d__1 = a[i__2].r;
#line 465 "chetf2_rook.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 465 "chetf2_rook.f"
		}
#line 467 "chetf2_rook.f"
	    }

/*           Update the leading submatrix */

#line 471 "chetf2_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

#line 479 "chetf2_rook.f"
		if (k > 1) {

/*                 Perform a rank-1 update of A(1:k-1,1:k-1) and */
/*                 store U(k) in column k */

#line 484 "chetf2_rook.f"
		    i__1 = k + k * a_dim1;
#line 484 "chetf2_rook.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) >= sfmin) {

/*                    Perform a rank-1 update of A(1:k-1,1:k-1) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*1/D(k)*W(k)**T */

#line 490 "chetf2_rook.f"
			i__1 = k + k * a_dim1;
#line 490 "chetf2_rook.f"
			d11 = 1. / a[i__1].r;
#line 491 "chetf2_rook.f"
			i__1 = k - 1;
#line 491 "chetf2_rook.f"
			d__1 = -d11;
#line 491 "chetf2_rook.f"
			cher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);

/*                    Store U(k) in column k */

#line 495 "chetf2_rook.f"
			i__1 = k - 1;
#line 495 "chetf2_rook.f"
			csscal_(&i__1, &d11, &a[k * a_dim1 + 1], &c__1);
#line 496 "chetf2_rook.f"
		    } else {

/*                    Store L(k) in column K */

#line 500 "chetf2_rook.f"
			i__1 = k + k * a_dim1;
#line 500 "chetf2_rook.f"
			d11 = a[i__1].r;
#line 501 "chetf2_rook.f"
			i__1 = k - 1;
#line 501 "chetf2_rook.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 502 "chetf2_rook.f"
			    i__2 = ii + k * a_dim1;
#line 502 "chetf2_rook.f"
			    i__3 = ii + k * a_dim1;
#line 502 "chetf2_rook.f"
			    z__1.r = a[i__3].r / d11, z__1.i = a[i__3].i / 
				    d11;
#line 502 "chetf2_rook.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 503 "chetf2_rook.f"
/* L16: */
#line 503 "chetf2_rook.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - U(k)*D(k)*U(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 510 "chetf2_rook.f"
			i__1 = k - 1;
#line 510 "chetf2_rook.f"
			d__1 = -d11;
#line 510 "chetf2_rook.f"
			cher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &
				a[a_offset], lda, (ftnlen)1);
#line 511 "chetf2_rook.f"
		    }
#line 512 "chetf2_rook.f"
		}

#line 514 "chetf2_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 530 "chetf2_rook.f"
		if (k > 2) {
/*                 D = |A12| */
#line 532 "chetf2_rook.f"
		    i__1 = k - 1 + k * a_dim1;
#line 532 "chetf2_rook.f"
		    d__1 = a[i__1].r;
#line 532 "chetf2_rook.f"
		    d__2 = d_imag(&a[k - 1 + k * a_dim1]);
#line 532 "chetf2_rook.f"
		    d__ = slapy2_(&d__1, &d__2);
#line 534 "chetf2_rook.f"
		    i__1 = k + k * a_dim1;
#line 534 "chetf2_rook.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 534 "chetf2_rook.f"
		    d11 = z__1.r;
#line 535 "chetf2_rook.f"
		    i__1 = k - 1 + (k - 1) * a_dim1;
#line 535 "chetf2_rook.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 535 "chetf2_rook.f"
		    d22 = z__1.r;
#line 536 "chetf2_rook.f"
		    i__1 = k - 1 + k * a_dim1;
#line 536 "chetf2_rook.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 536 "chetf2_rook.f"
		    d12.r = z__1.r, d12.i = z__1.i;
#line 537 "chetf2_rook.f"
		    tt = 1. / (d11 * d22 - 1.);

#line 539 "chetf2_rook.f"
		    for (j = k - 2; j >= 1; --j) {

/*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */

#line 543 "chetf2_rook.f"
			i__1 = j + (k - 1) * a_dim1;
#line 543 "chetf2_rook.f"
			z__3.r = d11 * a[i__1].r, z__3.i = d11 * a[i__1].i;
#line 543 "chetf2_rook.f"
			d_cnjg(&z__5, &d12);
#line 543 "chetf2_rook.f"
			i__2 = j + k * a_dim1;
#line 543 "chetf2_rook.f"
			z__4.r = z__5.r * a[i__2].r - z__5.i * a[i__2].i, 
				z__4.i = z__5.r * a[i__2].i + z__5.i * a[i__2]
				.r;
#line 543 "chetf2_rook.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 543 "chetf2_rook.f"
			z__1.r = tt * z__2.r, z__1.i = tt * z__2.i;
#line 543 "chetf2_rook.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 545 "chetf2_rook.f"
			i__1 = j + k * a_dim1;
#line 545 "chetf2_rook.f"
			z__3.r = d22 * a[i__1].r, z__3.i = d22 * a[i__1].i;
#line 545 "chetf2_rook.f"
			i__2 = j + (k - 1) * a_dim1;
#line 545 "chetf2_rook.f"
			z__4.r = d12.r * a[i__2].r - d12.i * a[i__2].i, 
				z__4.i = d12.r * a[i__2].i + d12.i * a[i__2]
				.r;
#line 545 "chetf2_rook.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 545 "chetf2_rook.f"
			z__1.r = tt * z__2.r, z__1.i = tt * z__2.i;
#line 545 "chetf2_rook.f"
			wk.r = z__1.r, wk.i = z__1.i;

/*                    Perform a rank-2 update of A(1:k-2,1:k-2) */

#line 549 "chetf2_rook.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 550 "chetf2_rook.f"
			    i__1 = i__ + j * a_dim1;
#line 550 "chetf2_rook.f"
			    i__2 = i__ + j * a_dim1;
#line 550 "chetf2_rook.f"
			    i__3 = i__ + k * a_dim1;
#line 550 "chetf2_rook.f"
			    z__4.r = a[i__3].r / d__, z__4.i = a[i__3].i / 
				    d__;
#line 550 "chetf2_rook.f"
			    d_cnjg(&z__5, &wk);
#line 550 "chetf2_rook.f"
			    z__3.r = z__4.r * z__5.r - z__4.i * z__5.i, 
				    z__3.i = z__4.r * z__5.i + z__4.i * 
				    z__5.r;
#line 550 "chetf2_rook.f"
			    z__2.r = a[i__2].r - z__3.r, z__2.i = a[i__2].i - 
				    z__3.i;
#line 550 "chetf2_rook.f"
			    i__4 = i__ + (k - 1) * a_dim1;
#line 550 "chetf2_rook.f"
			    z__7.r = a[i__4].r / d__, z__7.i = a[i__4].i / 
				    d__;
#line 550 "chetf2_rook.f"
			    d_cnjg(&z__8, &wkm1);
#line 550 "chetf2_rook.f"
			    z__6.r = z__7.r * z__8.r - z__7.i * z__8.i, 
				    z__6.i = z__7.r * z__8.i + z__7.i * 
				    z__8.r;
#line 550 "chetf2_rook.f"
			    z__1.r = z__2.r - z__6.r, z__1.i = z__2.i - 
				    z__6.i;
#line 550 "chetf2_rook.f"
			    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 553 "chetf2_rook.f"
/* L20: */
#line 553 "chetf2_rook.f"
			}

/*                    Store U(k) and U(k-1) in cols k and k-1 for row J */

#line 557 "chetf2_rook.f"
			i__1 = j + k * a_dim1;
#line 557 "chetf2_rook.f"
			z__1.r = wk.r / d__, z__1.i = wk.i / d__;
#line 557 "chetf2_rook.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 558 "chetf2_rook.f"
			i__1 = j + (k - 1) * a_dim1;
#line 558 "chetf2_rook.f"
			z__1.r = wkm1.r / d__, z__1.i = wkm1.i / d__;
#line 558 "chetf2_rook.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*                    (*) Make sure that diagonal element of pivot is real */
#line 560 "chetf2_rook.f"
			i__1 = j + j * a_dim1;
#line 560 "chetf2_rook.f"
			i__2 = j + j * a_dim1;
#line 560 "chetf2_rook.f"
			d__1 = a[i__2].r;
#line 560 "chetf2_rook.f"
			z__1.r = d__1, z__1.i = 0.;
#line 560 "chetf2_rook.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 562 "chetf2_rook.f"
/* L30: */
#line 562 "chetf2_rook.f"
		    }

#line 564 "chetf2_rook.f"
		}

#line 566 "chetf2_rook.f"
	    }

#line 568 "chetf2_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 572 "chetf2_rook.f"
	if (kstep == 1) {
#line 573 "chetf2_rook.f"
	    ipiv[k] = kp;
#line 574 "chetf2_rook.f"
	} else {
#line 575 "chetf2_rook.f"
	    ipiv[k] = -p;
#line 576 "chetf2_rook.f"
	    ipiv[k - 1] = -kp;
#line 577 "chetf2_rook.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 581 "chetf2_rook.f"
	k -= kstep;
#line 582 "chetf2_rook.f"
	goto L10;

#line 584 "chetf2_rook.f"
    } else {

/*        Factorize A as L*D*L**H using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 591 "chetf2_rook.f"
	k = 1;
#line 592 "chetf2_rook.f"
L40:

/*        If K > N, exit from loop */

#line 596 "chetf2_rook.f"
	if (k > *n) {
#line 596 "chetf2_rook.f"
	    goto L70;
#line 596 "chetf2_rook.f"
	}
#line 598 "chetf2_rook.f"
	kstep = 1;
#line 599 "chetf2_rook.f"
	p = k;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 604 "chetf2_rook.f"
	i__1 = k + k * a_dim1;
#line 604 "chetf2_rook.f"
	absakk = (d__1 = a[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 610 "chetf2_rook.f"
	if (k < *n) {
#line 611 "chetf2_rook.f"
	    i__1 = *n - k;
#line 611 "chetf2_rook.f"
	    imax = k + icamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 612 "chetf2_rook.f"
	    i__1 = imax + k * a_dim1;
#line 612 "chetf2_rook.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 613 "chetf2_rook.f"
	} else {
#line 614 "chetf2_rook.f"
	    colmax = 0.;
#line 615 "chetf2_rook.f"
	}

#line 617 "chetf2_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 621 "chetf2_rook.f"
	    if (*info == 0) {
#line 621 "chetf2_rook.f"
		*info = k;
#line 621 "chetf2_rook.f"
	    }
#line 623 "chetf2_rook.f"
	    kp = k;
#line 624 "chetf2_rook.f"
	    i__1 = k + k * a_dim1;
#line 624 "chetf2_rook.f"
	    i__2 = k + k * a_dim1;
#line 624 "chetf2_rook.f"
	    d__1 = a[i__2].r;
#line 624 "chetf2_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 625 "chetf2_rook.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 635 "chetf2_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 639 "chetf2_rook.f"
		kp = k;

#line 641 "chetf2_rook.f"
	    } else {

#line 643 "chetf2_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 647 "chetf2_rook.f"
L42:

/*                 BEGIN pivot search loop body */


/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 656 "chetf2_rook.f"
		if (imax != k) {
#line 657 "chetf2_rook.f"
		    i__1 = imax - k;
#line 657 "chetf2_rook.f"
		    jmax = k - 1 + icamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 658 "chetf2_rook.f"
		    i__1 = imax + jmax * a_dim1;
#line 658 "chetf2_rook.f"
		    rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    a[imax + jmax * a_dim1]), abs(d__2));
#line 659 "chetf2_rook.f"
		} else {
#line 660 "chetf2_rook.f"
		    rowmax = 0.;
#line 661 "chetf2_rook.f"
		}

#line 663 "chetf2_rook.f"
		if (imax < *n) {
#line 664 "chetf2_rook.f"
		    i__1 = *n - imax;
#line 664 "chetf2_rook.f"
		    itemp = imax + icamax_(&i__1, &a[imax + 1 + imax * a_dim1]
			    , &c__1);
#line 666 "chetf2_rook.f"
		    i__1 = itemp + imax * a_dim1;
#line 666 "chetf2_rook.f"
		    stemp = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    itemp + imax * a_dim1]), abs(d__2));
#line 667 "chetf2_rook.f"
		    if (stemp > rowmax) {
#line 668 "chetf2_rook.f"
			rowmax = stemp;
#line 669 "chetf2_rook.f"
			jmax = itemp;
#line 670 "chetf2_rook.f"
		    }
#line 671 "chetf2_rook.f"
		}

/*                 Case(2) */
/*                 Equivalent to testing for */
/*                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 678 "chetf2_rook.f"
		i__1 = imax + imax * a_dim1;
#line 678 "chetf2_rook.f"
		if (! ((d__1 = a[i__1].r, abs(d__1)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 684 "chetf2_rook.f"
		    kp = imax;
#line 685 "chetf2_rook.f"
		    done = TRUE_;

/*                 Case(3) */
/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 691 "chetf2_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 697 "chetf2_rook.f"
		    kp = imax;
#line 698 "chetf2_rook.f"
		    kstep = 2;
#line 699 "chetf2_rook.f"
		    done = TRUE_;

/*                 Case(4) */
#line 702 "chetf2_rook.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 706 "chetf2_rook.f"
		    p = imax;
#line 707 "chetf2_rook.f"
		    colmax = rowmax;
#line 708 "chetf2_rook.f"
		    imax = jmax;
#line 709 "chetf2_rook.f"
		}


/*                 END pivot search loop body */

#line 714 "chetf2_rook.f"
		if (! done) {
#line 714 "chetf2_rook.f"
		    goto L42;
#line 714 "chetf2_rook.f"
		}

#line 716 "chetf2_rook.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 724 "chetf2_rook.f"
	    kk = k + kstep - 1;

/*           For only a 2x2 pivot, interchange rows and columns K and P */
/*           in the trailing submatrix A(k:n,k:n) */

#line 729 "chetf2_rook.f"
	    if (kstep == 2 && p != k) {
/*              (1) Swap columnar parts */
#line 731 "chetf2_rook.f"
		if (p < *n) {
#line 731 "chetf2_rook.f"
		    i__1 = *n - p;
#line 731 "chetf2_rook.f"
		    cswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p 
			    * a_dim1], &c__1);
#line 731 "chetf2_rook.f"
		}
/*              (2) Swap and conjugate middle parts */
#line 734 "chetf2_rook.f"
		i__1 = p - 1;
#line 734 "chetf2_rook.f"
		for (j = k + 1; j <= i__1; ++j) {
#line 735 "chetf2_rook.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 735 "chetf2_rook.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 736 "chetf2_rook.f"
		    i__2 = j + k * a_dim1;
#line 736 "chetf2_rook.f"
		    d_cnjg(&z__1, &a[p + j * a_dim1]);
#line 736 "chetf2_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 737 "chetf2_rook.f"
		    i__2 = p + j * a_dim1;
#line 737 "chetf2_rook.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 738 "chetf2_rook.f"
/* L44: */
#line 738 "chetf2_rook.f"
		}
/*              (3) Swap and conjugate corner elements at row-col interserction */
#line 740 "chetf2_rook.f"
		i__1 = p + k * a_dim1;
#line 740 "chetf2_rook.f"
		d_cnjg(&z__1, &a[p + k * a_dim1]);
#line 740 "chetf2_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*              (4) Swap diagonal elements at row-col intersection */
#line 742 "chetf2_rook.f"
		i__1 = k + k * a_dim1;
#line 742 "chetf2_rook.f"
		r1 = a[i__1].r;
#line 743 "chetf2_rook.f"
		i__1 = k + k * a_dim1;
#line 743 "chetf2_rook.f"
		i__2 = p + p * a_dim1;
#line 743 "chetf2_rook.f"
		d__1 = a[i__2].r;
#line 743 "chetf2_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 744 "chetf2_rook.f"
		i__1 = p + p * a_dim1;
#line 744 "chetf2_rook.f"
		a[i__1].r = r1, a[i__1].i = 0.;
#line 745 "chetf2_rook.f"
	    }

/*           For both 1x1 and 2x2 pivots, interchange rows and */
/*           columns KK and KP in the trailing submatrix A(k:n,k:n) */

#line 750 "chetf2_rook.f"
	    if (kp != kk) {
/*              (1) Swap columnar parts */
#line 752 "chetf2_rook.f"
		if (kp < *n) {
#line 752 "chetf2_rook.f"
		    i__1 = *n - kp;
#line 752 "chetf2_rook.f"
		    cswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 752 "chetf2_rook.f"
		}
/*              (2) Swap and conjugate middle parts */
#line 755 "chetf2_rook.f"
		i__1 = kp - 1;
#line 755 "chetf2_rook.f"
		for (j = kk + 1; j <= i__1; ++j) {
#line 756 "chetf2_rook.f"
		    d_cnjg(&z__1, &a[j + kk * a_dim1]);
#line 756 "chetf2_rook.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 757 "chetf2_rook.f"
		    i__2 = j + kk * a_dim1;
#line 757 "chetf2_rook.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 757 "chetf2_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 758 "chetf2_rook.f"
		    i__2 = kp + j * a_dim1;
#line 758 "chetf2_rook.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 759 "chetf2_rook.f"
/* L45: */
#line 759 "chetf2_rook.f"
		}
/*              (3) Swap and conjugate corner elements at row-col interserction */
#line 761 "chetf2_rook.f"
		i__1 = kp + kk * a_dim1;
#line 761 "chetf2_rook.f"
		d_cnjg(&z__1, &a[kp + kk * a_dim1]);
#line 761 "chetf2_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
/*              (4) Swap diagonal elements at row-col intersection */
#line 763 "chetf2_rook.f"
		i__1 = kk + kk * a_dim1;
#line 763 "chetf2_rook.f"
		r1 = a[i__1].r;
#line 764 "chetf2_rook.f"
		i__1 = kk + kk * a_dim1;
#line 764 "chetf2_rook.f"
		i__2 = kp + kp * a_dim1;
#line 764 "chetf2_rook.f"
		d__1 = a[i__2].r;
#line 764 "chetf2_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 765 "chetf2_rook.f"
		i__1 = kp + kp * a_dim1;
#line 765 "chetf2_rook.f"
		a[i__1].r = r1, a[i__1].i = 0.;

#line 767 "chetf2_rook.f"
		if (kstep == 2) {
/*                 (*) Make sure that diagonal element of pivot is real */
#line 769 "chetf2_rook.f"
		    i__1 = k + k * a_dim1;
#line 769 "chetf2_rook.f"
		    i__2 = k + k * a_dim1;
#line 769 "chetf2_rook.f"
		    d__1 = a[i__2].r;
#line 769 "chetf2_rook.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
/*                 (5) Swap row elements */
#line 771 "chetf2_rook.f"
		    i__1 = k + 1 + k * a_dim1;
#line 771 "chetf2_rook.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 772 "chetf2_rook.f"
		    i__1 = k + 1 + k * a_dim1;
#line 772 "chetf2_rook.f"
		    i__2 = kp + k * a_dim1;
#line 772 "chetf2_rook.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 773 "chetf2_rook.f"
		    i__1 = kp + k * a_dim1;
#line 773 "chetf2_rook.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 774 "chetf2_rook.f"
		}
#line 775 "chetf2_rook.f"
	    } else {
/*              (*) Make sure that diagonal element of pivot is real */
#line 777 "chetf2_rook.f"
		i__1 = k + k * a_dim1;
#line 777 "chetf2_rook.f"
		i__2 = k + k * a_dim1;
#line 777 "chetf2_rook.f"
		d__1 = a[i__2].r;
#line 777 "chetf2_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 778 "chetf2_rook.f"
		if (kstep == 2) {
#line 778 "chetf2_rook.f"
		    i__1 = k + 1 + (k + 1) * a_dim1;
#line 778 "chetf2_rook.f"
		    i__2 = k + 1 + (k + 1) * a_dim1;
#line 778 "chetf2_rook.f"
		    d__1 = a[i__2].r;
#line 778 "chetf2_rook.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 778 "chetf2_rook.f"
		}
#line 780 "chetf2_rook.f"
	    }

/*           Update the trailing submatrix */

#line 784 "chetf2_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of A now holds */

/*              W(k) = L(k)*D(k), */

/*              where L(k) is the k-th column of L */

#line 792 "chetf2_rook.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) and */
/*                 store L(k) in column k */

/*                 Handle division by a small number */

#line 799 "chetf2_rook.f"
		    i__1 = k + k * a_dim1;
#line 799 "chetf2_rook.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) >= sfmin) {

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */

#line 805 "chetf2_rook.f"
			i__1 = k + k * a_dim1;
#line 805 "chetf2_rook.f"
			d11 = 1. / a[i__1].r;
#line 806 "chetf2_rook.f"
			i__1 = *n - k;
#line 806 "chetf2_rook.f"
			d__1 = -d11;
#line 806 "chetf2_rook.f"
			cher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);

/*                    Store L(k) in column k */

#line 811 "chetf2_rook.f"
			i__1 = *n - k;
#line 811 "chetf2_rook.f"
			csscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
#line 812 "chetf2_rook.f"
		    } else {

/*                    Store L(k) in column k */

#line 816 "chetf2_rook.f"
			i__1 = k + k * a_dim1;
#line 816 "chetf2_rook.f"
			d11 = a[i__1].r;
#line 817 "chetf2_rook.f"
			i__1 = *n;
#line 817 "chetf2_rook.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 818 "chetf2_rook.f"
			    i__2 = ii + k * a_dim1;
#line 818 "chetf2_rook.f"
			    i__3 = ii + k * a_dim1;
#line 818 "chetf2_rook.f"
			    z__1.r = a[i__3].r / d11, z__1.i = a[i__3].i / 
				    d11;
#line 818 "chetf2_rook.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 819 "chetf2_rook.f"
/* L46: */
#line 819 "chetf2_rook.f"
			}

/*                    Perform a rank-1 update of A(k+1:n,k+1:n) as */
/*                    A := A - L(k)*D(k)*L(k)**T */
/*                       = A - W(k)*(1/D(k))*W(k)**T */
/*                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */

#line 826 "chetf2_rook.f"
			i__1 = *n - k;
#line 826 "chetf2_rook.f"
			d__1 = -d11;
#line 826 "chetf2_rook.f"
			cher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &
				c__1, &a[k + 1 + (k + 1) * a_dim1], lda, (
				ftnlen)1);
#line 828 "chetf2_rook.f"
		    }
#line 829 "chetf2_rook.f"
		}

#line 831 "chetf2_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */


/*              Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T */
/*                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T */

/*              and store L(k) and L(k+1) in columns k and k+1 */

#line 848 "chetf2_rook.f"
		if (k < *n - 1) {
/*                 D = |A21| */
#line 850 "chetf2_rook.f"
		    i__1 = k + 1 + k * a_dim1;
#line 850 "chetf2_rook.f"
		    d__1 = a[i__1].r;
#line 850 "chetf2_rook.f"
		    d__2 = d_imag(&a[k + 1 + k * a_dim1]);
#line 850 "chetf2_rook.f"
		    d__ = slapy2_(&d__1, &d__2);
#line 852 "chetf2_rook.f"
		    i__1 = k + 1 + (k + 1) * a_dim1;
#line 852 "chetf2_rook.f"
		    d11 = a[i__1].r / d__;
#line 853 "chetf2_rook.f"
		    i__1 = k + k * a_dim1;
#line 853 "chetf2_rook.f"
		    d22 = a[i__1].r / d__;
#line 854 "chetf2_rook.f"
		    i__1 = k + 1 + k * a_dim1;
#line 854 "chetf2_rook.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 854 "chetf2_rook.f"
		    d21.r = z__1.r, d21.i = z__1.i;
#line 855 "chetf2_rook.f"
		    tt = 1. / (d11 * d22 - 1.);

#line 857 "chetf2_rook.f"
		    i__1 = *n;
#line 857 "chetf2_rook.f"
		    for (j = k + 2; j <= i__1; ++j) {

/*                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */

#line 861 "chetf2_rook.f"
			i__2 = j + k * a_dim1;
#line 861 "chetf2_rook.f"
			z__3.r = d11 * a[i__2].r, z__3.i = d11 * a[i__2].i;
#line 861 "chetf2_rook.f"
			i__3 = j + (k + 1) * a_dim1;
#line 861 "chetf2_rook.f"
			z__4.r = d21.r * a[i__3].r - d21.i * a[i__3].i, 
				z__4.i = d21.r * a[i__3].i + d21.i * a[i__3]
				.r;
#line 861 "chetf2_rook.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 861 "chetf2_rook.f"
			z__1.r = tt * z__2.r, z__1.i = tt * z__2.i;
#line 861 "chetf2_rook.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 862 "chetf2_rook.f"
			i__2 = j + (k + 1) * a_dim1;
#line 862 "chetf2_rook.f"
			z__3.r = d22 * a[i__2].r, z__3.i = d22 * a[i__2].i;
#line 862 "chetf2_rook.f"
			d_cnjg(&z__5, &d21);
#line 862 "chetf2_rook.f"
			i__3 = j + k * a_dim1;
#line 862 "chetf2_rook.f"
			z__4.r = z__5.r * a[i__3].r - z__5.i * a[i__3].i, 
				z__4.i = z__5.r * a[i__3].i + z__5.i * a[i__3]
				.r;
#line 862 "chetf2_rook.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 862 "chetf2_rook.f"
			z__1.r = tt * z__2.r, z__1.i = tt * z__2.i;
#line 862 "chetf2_rook.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;

/*                    Perform a rank-2 update of A(k+2:n,k+2:n) */

#line 867 "chetf2_rook.f"
			i__2 = *n;
#line 867 "chetf2_rook.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 868 "chetf2_rook.f"
			    i__3 = i__ + j * a_dim1;
#line 868 "chetf2_rook.f"
			    i__4 = i__ + j * a_dim1;
#line 868 "chetf2_rook.f"
			    i__5 = i__ + k * a_dim1;
#line 868 "chetf2_rook.f"
			    z__4.r = a[i__5].r / d__, z__4.i = a[i__5].i / 
				    d__;
#line 868 "chetf2_rook.f"
			    d_cnjg(&z__5, &wk);
#line 868 "chetf2_rook.f"
			    z__3.r = z__4.r * z__5.r - z__4.i * z__5.i, 
				    z__3.i = z__4.r * z__5.i + z__4.i * 
				    z__5.r;
#line 868 "chetf2_rook.f"
			    z__2.r = a[i__4].r - z__3.r, z__2.i = a[i__4].i - 
				    z__3.i;
#line 868 "chetf2_rook.f"
			    i__6 = i__ + (k + 1) * a_dim1;
#line 868 "chetf2_rook.f"
			    z__7.r = a[i__6].r / d__, z__7.i = a[i__6].i / 
				    d__;
#line 868 "chetf2_rook.f"
			    d_cnjg(&z__8, &wkp1);
#line 868 "chetf2_rook.f"
			    z__6.r = z__7.r * z__8.r - z__7.i * z__8.i, 
				    z__6.i = z__7.r * z__8.i + z__7.i * 
				    z__8.r;
#line 868 "chetf2_rook.f"
			    z__1.r = z__2.r - z__6.r, z__1.i = z__2.i - 
				    z__6.i;
#line 868 "chetf2_rook.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 871 "chetf2_rook.f"
/* L50: */
#line 871 "chetf2_rook.f"
			}

/*                    Store L(k) and L(k+1) in cols k and k+1 for row J */

#line 875 "chetf2_rook.f"
			i__2 = j + k * a_dim1;
#line 875 "chetf2_rook.f"
			z__1.r = wk.r / d__, z__1.i = wk.i / d__;
#line 875 "chetf2_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 876 "chetf2_rook.f"
			i__2 = j + (k + 1) * a_dim1;
#line 876 "chetf2_rook.f"
			z__1.r = wkp1.r / d__, z__1.i = wkp1.i / d__;
#line 876 "chetf2_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/*                    (*) Make sure that diagonal element of pivot is real */
#line 878 "chetf2_rook.f"
			i__2 = j + j * a_dim1;
#line 878 "chetf2_rook.f"
			i__3 = j + j * a_dim1;
#line 878 "chetf2_rook.f"
			d__1 = a[i__3].r;
#line 878 "chetf2_rook.f"
			z__1.r = d__1, z__1.i = 0.;
#line 878 "chetf2_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;

#line 880 "chetf2_rook.f"
/* L60: */
#line 880 "chetf2_rook.f"
		    }

#line 882 "chetf2_rook.f"
		}

#line 884 "chetf2_rook.f"
	    }

#line 886 "chetf2_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 890 "chetf2_rook.f"
	if (kstep == 1) {
#line 891 "chetf2_rook.f"
	    ipiv[k] = kp;
#line 892 "chetf2_rook.f"
	} else {
#line 893 "chetf2_rook.f"
	    ipiv[k] = -p;
#line 894 "chetf2_rook.f"
	    ipiv[k + 1] = -kp;
#line 895 "chetf2_rook.f"
	}

/*        Increase K and return to the start of the main loop */

#line 899 "chetf2_rook.f"
	k += kstep;
#line 900 "chetf2_rook.f"
	goto L40;

#line 902 "chetf2_rook.f"
    }

#line 904 "chetf2_rook.f"
L70:

#line 906 "chetf2_rook.f"
    return 0;

/*     End of CHETF2_ROOK */

} /* chetf2_rook__ */

