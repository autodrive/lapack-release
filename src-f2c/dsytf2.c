#line 1 "dsytf2.f"
/* dsytf2.f -- translated by f2c (version 20100827).
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

#line 1 "dsytf2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal piv
oting method (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTF2 computes the factorization of a real symmetric matrix A using */
/* > the Bunch-Kaufman diagonal pivoting method: */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >             interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) = IPIV(k-1) < 0, then rows and columns */
/* >             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* >             is a 2-by-2 diagonal block. */
/* > */
/* >          If UPLO = 'L': */
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >             interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) = IPIV(k+1) < 0, then rows and columns */
/* >             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1) */
/* >             is a 2-by-2 diagonal block. */
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

/* > \date December 2016 */

/* > \ingroup doubleSYcomputational */

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
/* >  09-29-06 - patch from */
/* >    Bobby Cheng, MathWorks */
/* > */
/* >    Replace l.204 and l.372 */
/* >         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN */
/* >    by */
/* >         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN */
/* > */
/* >  01-01-96 - Based on modifications by */
/* >    J. Lewis, Boeing Computer Services Company */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* >  1-96 - Based on modifications by J. Lewis, Boeing Computer Services */
/* >         Company */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int dsytf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal t, r1, d11, d12, d21, d22;
    static integer kk, kp;
    static doublereal wk, wkm1, wkp1;
    static integer imax, jmax;
    extern /* Subroutine */ int dsyr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer kstep;
    static logical upper;
    static doublereal absakk;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern logical disnan_(doublereal *);
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

#line 240 "dsytf2.f"
    /* Parameter adjustments */
#line 240 "dsytf2.f"
    a_dim1 = *lda;
#line 240 "dsytf2.f"
    a_offset = 1 + a_dim1;
#line 240 "dsytf2.f"
    a -= a_offset;
#line 240 "dsytf2.f"
    --ipiv;
#line 240 "dsytf2.f"

#line 240 "dsytf2.f"
    /* Function Body */
#line 240 "dsytf2.f"
    *info = 0;
#line 241 "dsytf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 242 "dsytf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 243 "dsytf2.f"
	*info = -1;
#line 244 "dsytf2.f"
    } else if (*n < 0) {
#line 245 "dsytf2.f"
	*info = -2;
#line 246 "dsytf2.f"
    } else if (*lda < max(1,*n)) {
#line 247 "dsytf2.f"
	*info = -4;
#line 248 "dsytf2.f"
    }
#line 249 "dsytf2.f"
    if (*info != 0) {
#line 250 "dsytf2.f"
	i__1 = -(*info);
#line 250 "dsytf2.f"
	xerbla_("DSYTF2", &i__1, (ftnlen)6);
#line 251 "dsytf2.f"
	return 0;
#line 252 "dsytf2.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 256 "dsytf2.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 258 "dsytf2.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 265 "dsytf2.f"
	k = *n;
#line 266 "dsytf2.f"
L10:

/*        If K < 1, exit from loop */

#line 270 "dsytf2.f"
	if (k < 1) {
#line 270 "dsytf2.f"
	    goto L70;
#line 270 "dsytf2.f"
	}
#line 272 "dsytf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 277 "dsytf2.f"
	absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 283 "dsytf2.f"
	if (k > 1) {
#line 284 "dsytf2.f"
	    i__1 = k - 1;
#line 284 "dsytf2.f"
	    imax = idamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 285 "dsytf2.f"
	    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
#line 286 "dsytf2.f"
	} else {
#line 287 "dsytf2.f"
	    colmax = 0.;
#line 288 "dsytf2.f"
	}

#line 290 "dsytf2.f"
	if (max(absakk,colmax) == 0. || disnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 295 "dsytf2.f"
	    if (*info == 0) {
#line 295 "dsytf2.f"
		*info = k;
#line 295 "dsytf2.f"
	    }
#line 297 "dsytf2.f"
	    kp = k;
#line 298 "dsytf2.f"
	} else {
#line 299 "dsytf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 303 "dsytf2.f"
		kp = k;
#line 304 "dsytf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 309 "dsytf2.f"
		i__1 = k - imax;
#line 309 "dsytf2.f"
		jmax = imax + idamax_(&i__1, &a[imax + (imax + 1) * a_dim1], 
			lda);
#line 310 "dsytf2.f"
		rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
#line 311 "dsytf2.f"
		if (imax > 1) {
#line 312 "dsytf2.f"
		    i__1 = imax - 1;
#line 312 "dsytf2.f"
		    jmax = idamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
#line 313 "dsytf2.f"
		    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], 
			    abs(d__1));
#line 313 "dsytf2.f"
		    rowmax = max(d__2,d__3);
#line 314 "dsytf2.f"
		}

#line 316 "dsytf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 320 "dsytf2.f"
		    kp = k;
#line 321 "dsytf2.f"
		} else if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) >= 
			alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 326 "dsytf2.f"
		    kp = imax;
#line 327 "dsytf2.f"
		} else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 332 "dsytf2.f"
		    kp = imax;
#line 333 "dsytf2.f"
		    kstep = 2;
#line 334 "dsytf2.f"
		}
#line 335 "dsytf2.f"
	    }

#line 337 "dsytf2.f"
	    kk = k - kstep + 1;
#line 338 "dsytf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 343 "dsytf2.f"
		i__1 = kp - 1;
#line 343 "dsytf2.f"
		dswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1],
			 &c__1);
#line 344 "dsytf2.f"
		i__1 = kk - kp - 1;
#line 344 "dsytf2.f"
		dswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 346 "dsytf2.f"
		t = a[kk + kk * a_dim1];
#line 347 "dsytf2.f"
		a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
#line 348 "dsytf2.f"
		a[kp + kp * a_dim1] = t;
#line 349 "dsytf2.f"
		if (kstep == 2) {
#line 350 "dsytf2.f"
		    t = a[k - 1 + k * a_dim1];
#line 351 "dsytf2.f"
		    a[k - 1 + k * a_dim1] = a[kp + k * a_dim1];
#line 352 "dsytf2.f"
		    a[kp + k * a_dim1] = t;
#line 353 "dsytf2.f"
		}
#line 354 "dsytf2.f"
	    }

/*           Update the leading submatrix */

#line 358 "dsytf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */

#line 370 "dsytf2.f"
		r1 = 1. / a[k + k * a_dim1];
#line 371 "dsytf2.f"
		i__1 = k - 1;
#line 371 "dsytf2.f"
		d__1 = -r1;
#line 371 "dsytf2.f"
		dsyr_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &a[
			a_offset], lda, (ftnlen)1);

/*              Store U(k) in column k */

#line 375 "dsytf2.f"
		i__1 = k - 1;
#line 375 "dsytf2.f"
		dscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 376 "dsytf2.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */

#line 390 "dsytf2.f"
		if (k > 2) {

#line 392 "dsytf2.f"
		    d12 = a[k - 1 + k * a_dim1];
#line 393 "dsytf2.f"
		    d22 = a[k - 1 + (k - 1) * a_dim1] / d12;
#line 394 "dsytf2.f"
		    d11 = a[k + k * a_dim1] / d12;
#line 395 "dsytf2.f"
		    t = 1. / (d11 * d22 - 1.);
#line 396 "dsytf2.f"
		    d12 = t / d12;

#line 398 "dsytf2.f"
		    for (j = k - 2; j >= 1; --j) {
#line 399 "dsytf2.f"
			wkm1 = d12 * (d11 * a[j + (k - 1) * a_dim1] - a[j + k 
				* a_dim1]);
#line 400 "dsytf2.f"
			wk = d12 * (d22 * a[j + k * a_dim1] - a[j + (k - 1) * 
				a_dim1]);
#line 401 "dsytf2.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 402 "dsytf2.f"
			    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ 
				    + k * a_dim1] * wk - a[i__ + (k - 1) * 
				    a_dim1] * wkm1;
#line 404 "dsytf2.f"
/* L20: */
#line 404 "dsytf2.f"
			}
#line 405 "dsytf2.f"
			a[j + k * a_dim1] = wk;
#line 406 "dsytf2.f"
			a[j + (k - 1) * a_dim1] = wkm1;
#line 407 "dsytf2.f"
/* L30: */
#line 407 "dsytf2.f"
		    }

#line 409 "dsytf2.f"
		}

#line 411 "dsytf2.f"
	    }
#line 412 "dsytf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 416 "dsytf2.f"
	if (kstep == 1) {
#line 417 "dsytf2.f"
	    ipiv[k] = kp;
#line 418 "dsytf2.f"
	} else {
#line 419 "dsytf2.f"
	    ipiv[k] = -kp;
#line 420 "dsytf2.f"
	    ipiv[k - 1] = -kp;
#line 421 "dsytf2.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 425 "dsytf2.f"
	k -= kstep;
#line 426 "dsytf2.f"
	goto L10;

#line 428 "dsytf2.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 435 "dsytf2.f"
	k = 1;
#line 436 "dsytf2.f"
L40:

/*        If K > N, exit from loop */

#line 440 "dsytf2.f"
	if (k > *n) {
#line 440 "dsytf2.f"
	    goto L70;
#line 440 "dsytf2.f"
	}
#line 442 "dsytf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 447 "dsytf2.f"
	absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 453 "dsytf2.f"
	if (k < *n) {
#line 454 "dsytf2.f"
	    i__1 = *n - k;
#line 454 "dsytf2.f"
	    imax = k + idamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 455 "dsytf2.f"
	    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
#line 456 "dsytf2.f"
	} else {
#line 457 "dsytf2.f"
	    colmax = 0.;
#line 458 "dsytf2.f"
	}

#line 460 "dsytf2.f"
	if (max(absakk,colmax) == 0. || disnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 465 "dsytf2.f"
	    if (*info == 0) {
#line 465 "dsytf2.f"
		*info = k;
#line 465 "dsytf2.f"
	    }
#line 467 "dsytf2.f"
	    kp = k;
#line 468 "dsytf2.f"
	} else {
#line 469 "dsytf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 473 "dsytf2.f"
		kp = k;
#line 474 "dsytf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 479 "dsytf2.f"
		i__1 = imax - k;
#line 479 "dsytf2.f"
		jmax = k - 1 + idamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 480 "dsytf2.f"
		rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
#line 481 "dsytf2.f"
		if (imax < *n) {
#line 482 "dsytf2.f"
		    i__1 = *n - imax;
#line 482 "dsytf2.f"
		    jmax = imax + idamax_(&i__1, &a[imax + 1 + imax * a_dim1],
			     &c__1);
/* Computing MAX */
#line 483 "dsytf2.f"
		    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], 
			    abs(d__1));
#line 483 "dsytf2.f"
		    rowmax = max(d__2,d__3);
#line 484 "dsytf2.f"
		}

#line 486 "dsytf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 490 "dsytf2.f"
		    kp = k;
#line 491 "dsytf2.f"
		} else if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) >= 
			alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 496 "dsytf2.f"
		    kp = imax;
#line 497 "dsytf2.f"
		} else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 502 "dsytf2.f"
		    kp = imax;
#line 503 "dsytf2.f"
		    kstep = 2;
#line 504 "dsytf2.f"
		}
#line 505 "dsytf2.f"
	    }

#line 507 "dsytf2.f"
	    kk = k + kstep - 1;
#line 508 "dsytf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 513 "dsytf2.f"
		if (kp < *n) {
#line 513 "dsytf2.f"
		    i__1 = *n - kp;
#line 513 "dsytf2.f"
		    dswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 513 "dsytf2.f"
		}
#line 515 "dsytf2.f"
		i__1 = kp - kk - 1;
#line 515 "dsytf2.f"
		dswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 517 "dsytf2.f"
		t = a[kk + kk * a_dim1];
#line 518 "dsytf2.f"
		a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
#line 519 "dsytf2.f"
		a[kp + kp * a_dim1] = t;
#line 520 "dsytf2.f"
		if (kstep == 2) {
#line 521 "dsytf2.f"
		    t = a[k + 1 + k * a_dim1];
#line 522 "dsytf2.f"
		    a[k + 1 + k * a_dim1] = a[kp + k * a_dim1];
#line 523 "dsytf2.f"
		    a[kp + k * a_dim1] = t;
#line 524 "dsytf2.f"
		}
#line 525 "dsytf2.f"
	    }

/*           Update the trailing submatrix */

#line 529 "dsytf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 537 "dsytf2.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */

#line 543 "dsytf2.f"
		    d11 = 1. / a[k + k * a_dim1];
#line 544 "dsytf2.f"
		    i__1 = *n - k;
#line 544 "dsytf2.f"
		    d__1 = -d11;
#line 544 "dsytf2.f"
		    dsyr_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &c__1, &
			    a[k + 1 + (k + 1) * a_dim1], lda, (ftnlen)1);

/*                 Store L(k) in column K */

#line 549 "dsytf2.f"
		    i__1 = *n - k;
#line 549 "dsytf2.f"
		    dscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
#line 550 "dsytf2.f"
		}
#line 551 "dsytf2.f"
	    } else {

/*              2-by-2 pivot block D(k) */

#line 555 "dsytf2.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))**T */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 564 "dsytf2.f"
		    d21 = a[k + 1 + k * a_dim1];
#line 565 "dsytf2.f"
		    d11 = a[k + 1 + (k + 1) * a_dim1] / d21;
#line 566 "dsytf2.f"
		    d22 = a[k + k * a_dim1] / d21;
#line 567 "dsytf2.f"
		    t = 1. / (d11 * d22 - 1.);
#line 568 "dsytf2.f"
		    d21 = t / d21;

#line 570 "dsytf2.f"
		    i__1 = *n;
#line 570 "dsytf2.f"
		    for (j = k + 2; j <= i__1; ++j) {

#line 572 "dsytf2.f"
			wk = d21 * (d11 * a[j + k * a_dim1] - a[j + (k + 1) * 
				a_dim1]);
#line 573 "dsytf2.f"
			wkp1 = d21 * (d22 * a[j + (k + 1) * a_dim1] - a[j + k 
				* a_dim1]);

#line 575 "dsytf2.f"
			i__2 = *n;
#line 575 "dsytf2.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 576 "dsytf2.f"
			    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ 
				    + k * a_dim1] * wk - a[i__ + (k + 1) * 
				    a_dim1] * wkp1;
#line 578 "dsytf2.f"
/* L50: */
#line 578 "dsytf2.f"
			}

#line 580 "dsytf2.f"
			a[j + k * a_dim1] = wk;
#line 581 "dsytf2.f"
			a[j + (k + 1) * a_dim1] = wkp1;

#line 583 "dsytf2.f"
/* L60: */
#line 583 "dsytf2.f"
		    }
#line 584 "dsytf2.f"
		}
#line 585 "dsytf2.f"
	    }
#line 586 "dsytf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 590 "dsytf2.f"
	if (kstep == 1) {
#line 591 "dsytf2.f"
	    ipiv[k] = kp;
#line 592 "dsytf2.f"
	} else {
#line 593 "dsytf2.f"
	    ipiv[k] = -kp;
#line 594 "dsytf2.f"
	    ipiv[k + 1] = -kp;
#line 595 "dsytf2.f"
	}

/*        Increase K and return to the start of the main loop */

#line 599 "dsytf2.f"
	k += kstep;
#line 600 "dsytf2.f"
	goto L40;

#line 602 "dsytf2.f"
    }

#line 604 "dsytf2.f"
L70:

#line 606 "dsytf2.f"
    return 0;

/*     End of DSYTF2 */

} /* dsytf2_ */

