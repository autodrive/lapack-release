#line 1 "ssytf2.f"
/* ssytf2.f -- translated by f2c (version 20100827).
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

#line 1 "ssytf2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal piv
oting method (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTF2( UPLO, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTF2 computes the factorization of a real symmetric matrix A using */
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
/* >          A is REAL array, dimension (LDA,N) */
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

/* > \ingroup realSYcomputational */

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
/* >         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN */
/* > */
/* >  01-01-96 - Based on modifications by */
/* >    J. Lewis, Boeing Computer Services Company */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* >  1-96 - Based on modifications by J. Lewis, Boeing Computer Services */
/* >         Company */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int ssytf2_(char *uplo, integer *n, doublereal *a, integer *
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
    extern /* Subroutine */ int ssyr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal absakk;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal colmax;
    extern logical sisnan_(doublereal *);
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 241 "ssytf2.f"
    /* Parameter adjustments */
#line 241 "ssytf2.f"
    a_dim1 = *lda;
#line 241 "ssytf2.f"
    a_offset = 1 + a_dim1;
#line 241 "ssytf2.f"
    a -= a_offset;
#line 241 "ssytf2.f"
    --ipiv;
#line 241 "ssytf2.f"

#line 241 "ssytf2.f"
    /* Function Body */
#line 241 "ssytf2.f"
    *info = 0;
#line 242 "ssytf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 243 "ssytf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 244 "ssytf2.f"
	*info = -1;
#line 245 "ssytf2.f"
    } else if (*n < 0) {
#line 246 "ssytf2.f"
	*info = -2;
#line 247 "ssytf2.f"
    } else if (*lda < max(1,*n)) {
#line 248 "ssytf2.f"
	*info = -4;
#line 249 "ssytf2.f"
    }
#line 250 "ssytf2.f"
    if (*info != 0) {
#line 251 "ssytf2.f"
	i__1 = -(*info);
#line 251 "ssytf2.f"
	xerbla_("SSYTF2", &i__1, (ftnlen)6);
#line 252 "ssytf2.f"
	return 0;
#line 253 "ssytf2.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 257 "ssytf2.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 259 "ssytf2.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 266 "ssytf2.f"
	k = *n;
#line 267 "ssytf2.f"
L10:

/*        If K < 1, exit from loop */

#line 271 "ssytf2.f"
	if (k < 1) {
#line 271 "ssytf2.f"
	    goto L70;
#line 271 "ssytf2.f"
	}
#line 273 "ssytf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 278 "ssytf2.f"
	absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 284 "ssytf2.f"
	if (k > 1) {
#line 285 "ssytf2.f"
	    i__1 = k - 1;
#line 285 "ssytf2.f"
	    imax = isamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 286 "ssytf2.f"
	    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
#line 287 "ssytf2.f"
	} else {
#line 288 "ssytf2.f"
	    colmax = 0.;
#line 289 "ssytf2.f"
	}

#line 291 "ssytf2.f"
	if (max(absakk,colmax) == 0. || sisnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 296 "ssytf2.f"
	    if (*info == 0) {
#line 296 "ssytf2.f"
		*info = k;
#line 296 "ssytf2.f"
	    }
#line 298 "ssytf2.f"
	    kp = k;
#line 299 "ssytf2.f"
	} else {
#line 300 "ssytf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 304 "ssytf2.f"
		kp = k;
#line 305 "ssytf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 310 "ssytf2.f"
		i__1 = k - imax;
#line 310 "ssytf2.f"
		jmax = imax + isamax_(&i__1, &a[imax + (imax + 1) * a_dim1], 
			lda);
#line 311 "ssytf2.f"
		rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
#line 312 "ssytf2.f"
		if (imax > 1) {
#line 313 "ssytf2.f"
		    i__1 = imax - 1;
#line 313 "ssytf2.f"
		    jmax = isamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
#line 314 "ssytf2.f"
		    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], 
			    abs(d__1));
#line 314 "ssytf2.f"
		    rowmax = max(d__2,d__3);
#line 315 "ssytf2.f"
		}

#line 317 "ssytf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 321 "ssytf2.f"
		    kp = k;
#line 322 "ssytf2.f"
		} else if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) >= 
			alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 327 "ssytf2.f"
		    kp = imax;
#line 328 "ssytf2.f"
		} else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 333 "ssytf2.f"
		    kp = imax;
#line 334 "ssytf2.f"
		    kstep = 2;
#line 335 "ssytf2.f"
		}
#line 336 "ssytf2.f"
	    }

#line 338 "ssytf2.f"
	    kk = k - kstep + 1;
#line 339 "ssytf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 344 "ssytf2.f"
		i__1 = kp - 1;
#line 344 "ssytf2.f"
		sswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1],
			 &c__1);
#line 345 "ssytf2.f"
		i__1 = kk - kp - 1;
#line 345 "ssytf2.f"
		sswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 347 "ssytf2.f"
		t = a[kk + kk * a_dim1];
#line 348 "ssytf2.f"
		a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
#line 349 "ssytf2.f"
		a[kp + kp * a_dim1] = t;
#line 350 "ssytf2.f"
		if (kstep == 2) {
#line 351 "ssytf2.f"
		    t = a[k - 1 + k * a_dim1];
#line 352 "ssytf2.f"
		    a[k - 1 + k * a_dim1] = a[kp + k * a_dim1];
#line 353 "ssytf2.f"
		    a[kp + k * a_dim1] = t;
#line 354 "ssytf2.f"
		}
#line 355 "ssytf2.f"
	    }

/*           Update the leading submatrix */

#line 359 "ssytf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */

#line 371 "ssytf2.f"
		r1 = 1. / a[k + k * a_dim1];
#line 372 "ssytf2.f"
		i__1 = k - 1;
#line 372 "ssytf2.f"
		d__1 = -r1;
#line 372 "ssytf2.f"
		ssyr_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &a[
			a_offset], lda, (ftnlen)1);

/*              Store U(k) in column k */

#line 376 "ssytf2.f"
		i__1 = k - 1;
#line 376 "ssytf2.f"
		sscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 377 "ssytf2.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */

#line 391 "ssytf2.f"
		if (k > 2) {

#line 393 "ssytf2.f"
		    d12 = a[k - 1 + k * a_dim1];
#line 394 "ssytf2.f"
		    d22 = a[k - 1 + (k - 1) * a_dim1] / d12;
#line 395 "ssytf2.f"
		    d11 = a[k + k * a_dim1] / d12;
#line 396 "ssytf2.f"
		    t = 1. / (d11 * d22 - 1.);
#line 397 "ssytf2.f"
		    d12 = t / d12;

#line 399 "ssytf2.f"
		    for (j = k - 2; j >= 1; --j) {
#line 400 "ssytf2.f"
			wkm1 = d12 * (d11 * a[j + (k - 1) * a_dim1] - a[j + k 
				* a_dim1]);
#line 401 "ssytf2.f"
			wk = d12 * (d22 * a[j + k * a_dim1] - a[j + (k - 1) * 
				a_dim1]);
#line 402 "ssytf2.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 403 "ssytf2.f"
			    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ 
				    + k * a_dim1] * wk - a[i__ + (k - 1) * 
				    a_dim1] * wkm1;
#line 405 "ssytf2.f"
/* L20: */
#line 405 "ssytf2.f"
			}
#line 406 "ssytf2.f"
			a[j + k * a_dim1] = wk;
#line 407 "ssytf2.f"
			a[j + (k - 1) * a_dim1] = wkm1;
#line 408 "ssytf2.f"
/* L30: */
#line 408 "ssytf2.f"
		    }

#line 410 "ssytf2.f"
		}

#line 412 "ssytf2.f"
	    }
#line 413 "ssytf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 417 "ssytf2.f"
	if (kstep == 1) {
#line 418 "ssytf2.f"
	    ipiv[k] = kp;
#line 419 "ssytf2.f"
	} else {
#line 420 "ssytf2.f"
	    ipiv[k] = -kp;
#line 421 "ssytf2.f"
	    ipiv[k - 1] = -kp;
#line 422 "ssytf2.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 426 "ssytf2.f"
	k -= kstep;
#line 427 "ssytf2.f"
	goto L10;

#line 429 "ssytf2.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 436 "ssytf2.f"
	k = 1;
#line 437 "ssytf2.f"
L40:

/*        If K > N, exit from loop */

#line 441 "ssytf2.f"
	if (k > *n) {
#line 441 "ssytf2.f"
	    goto L70;
#line 441 "ssytf2.f"
	}
#line 443 "ssytf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 448 "ssytf2.f"
	absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 454 "ssytf2.f"
	if (k < *n) {
#line 455 "ssytf2.f"
	    i__1 = *n - k;
#line 455 "ssytf2.f"
	    imax = k + isamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 456 "ssytf2.f"
	    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
#line 457 "ssytf2.f"
	} else {
#line 458 "ssytf2.f"
	    colmax = 0.;
#line 459 "ssytf2.f"
	}

#line 461 "ssytf2.f"
	if (max(absakk,colmax) == 0. || sisnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 466 "ssytf2.f"
	    if (*info == 0) {
#line 466 "ssytf2.f"
		*info = k;
#line 466 "ssytf2.f"
	    }
#line 468 "ssytf2.f"
	    kp = k;
#line 469 "ssytf2.f"
	} else {
#line 470 "ssytf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 474 "ssytf2.f"
		kp = k;
#line 475 "ssytf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 480 "ssytf2.f"
		i__1 = imax - k;
#line 480 "ssytf2.f"
		jmax = k - 1 + isamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 481 "ssytf2.f"
		rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
#line 482 "ssytf2.f"
		if (imax < *n) {
#line 483 "ssytf2.f"
		    i__1 = *n - imax;
#line 483 "ssytf2.f"
		    jmax = imax + isamax_(&i__1, &a[imax + 1 + imax * a_dim1],
			     &c__1);
/* Computing MAX */
#line 484 "ssytf2.f"
		    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], 
			    abs(d__1));
#line 484 "ssytf2.f"
		    rowmax = max(d__2,d__3);
#line 485 "ssytf2.f"
		}

#line 487 "ssytf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 491 "ssytf2.f"
		    kp = k;
#line 492 "ssytf2.f"
		} else if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) >= 
			alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 497 "ssytf2.f"
		    kp = imax;
#line 498 "ssytf2.f"
		} else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 503 "ssytf2.f"
		    kp = imax;
#line 504 "ssytf2.f"
		    kstep = 2;
#line 505 "ssytf2.f"
		}
#line 506 "ssytf2.f"
	    }

#line 508 "ssytf2.f"
	    kk = k + kstep - 1;
#line 509 "ssytf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 514 "ssytf2.f"
		if (kp < *n) {
#line 514 "ssytf2.f"
		    i__1 = *n - kp;
#line 514 "ssytf2.f"
		    sswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 514 "ssytf2.f"
		}
#line 516 "ssytf2.f"
		i__1 = kp - kk - 1;
#line 516 "ssytf2.f"
		sswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 518 "ssytf2.f"
		t = a[kk + kk * a_dim1];
#line 519 "ssytf2.f"
		a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
#line 520 "ssytf2.f"
		a[kp + kp * a_dim1] = t;
#line 521 "ssytf2.f"
		if (kstep == 2) {
#line 522 "ssytf2.f"
		    t = a[k + 1 + k * a_dim1];
#line 523 "ssytf2.f"
		    a[k + 1 + k * a_dim1] = a[kp + k * a_dim1];
#line 524 "ssytf2.f"
		    a[kp + k * a_dim1] = t;
#line 525 "ssytf2.f"
		}
#line 526 "ssytf2.f"
	    }

/*           Update the trailing submatrix */

#line 530 "ssytf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 538 "ssytf2.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */

#line 544 "ssytf2.f"
		    d11 = 1. / a[k + k * a_dim1];
#line 545 "ssytf2.f"
		    i__1 = *n - k;
#line 545 "ssytf2.f"
		    d__1 = -d11;
#line 545 "ssytf2.f"
		    ssyr_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &c__1, &
			    a[k + 1 + (k + 1) * a_dim1], lda, (ftnlen)1);

/*                 Store L(k) in column K */

#line 550 "ssytf2.f"
		    i__1 = *n - k;
#line 550 "ssytf2.f"
		    sscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
#line 551 "ssytf2.f"
		}
#line 552 "ssytf2.f"
	    } else {

/*              2-by-2 pivot block D(k) */

#line 556 "ssytf2.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))**T */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 565 "ssytf2.f"
		    d21 = a[k + 1 + k * a_dim1];
#line 566 "ssytf2.f"
		    d11 = a[k + 1 + (k + 1) * a_dim1] / d21;
#line 567 "ssytf2.f"
		    d22 = a[k + k * a_dim1] / d21;
#line 568 "ssytf2.f"
		    t = 1. / (d11 * d22 - 1.);
#line 569 "ssytf2.f"
		    d21 = t / d21;

#line 571 "ssytf2.f"
		    i__1 = *n;
#line 571 "ssytf2.f"
		    for (j = k + 2; j <= i__1; ++j) {

#line 573 "ssytf2.f"
			wk = d21 * (d11 * a[j + k * a_dim1] - a[j + (k + 1) * 
				a_dim1]);
#line 574 "ssytf2.f"
			wkp1 = d21 * (d22 * a[j + (k + 1) * a_dim1] - a[j + k 
				* a_dim1]);

#line 576 "ssytf2.f"
			i__2 = *n;
#line 576 "ssytf2.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 577 "ssytf2.f"
			    a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ 
				    + k * a_dim1] * wk - a[i__ + (k + 1) * 
				    a_dim1] * wkp1;
#line 579 "ssytf2.f"
/* L50: */
#line 579 "ssytf2.f"
			}

#line 581 "ssytf2.f"
			a[j + k * a_dim1] = wk;
#line 582 "ssytf2.f"
			a[j + (k + 1) * a_dim1] = wkp1;

#line 584 "ssytf2.f"
/* L60: */
#line 584 "ssytf2.f"
		    }
#line 585 "ssytf2.f"
		}
#line 586 "ssytf2.f"
	    }
#line 587 "ssytf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 591 "ssytf2.f"
	if (kstep == 1) {
#line 592 "ssytf2.f"
	    ipiv[k] = kp;
#line 593 "ssytf2.f"
	} else {
#line 594 "ssytf2.f"
	    ipiv[k] = -kp;
#line 595 "ssytf2.f"
	    ipiv[k + 1] = -kp;
#line 596 "ssytf2.f"
	}

/*        Increase K and return to the start of the main loop */

#line 600 "ssytf2.f"
	k += kstep;
#line 601 "ssytf2.f"
	goto L40;

#line 603 "ssytf2.f"
    }

#line 605 "ssytf2.f"
L70:

#line 607 "ssytf2.f"
    return 0;

/*     End of SSYTF2 */

} /* ssytf2_ */

