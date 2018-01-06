#line 1 "dsptrf.f"
/* dsptrf.f -- translated by f2c (version 20100827).
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

#line 1 "dsptrf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DSPTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsptrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsptrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsptrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPTRF computes the factorization of a real symmetric matrix A stored */
/* > in packed format using the Bunch-Kaufman diagonal pivoting method: */
/* > */
/* >    A = U*D*U**T  or  A = L*D*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is symmetric and block diagonal with */
/* > 1-by-1 and 2-by-2 diagonal blocks. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the symmetric matrix */
/* >          A, packed columnwise in a linear array.  The j-th column of A */
/* >          is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* >          On exit, the block diagonal matrix D and the multipliers used */
/* >          to obtain the factor U or L, stored as a packed triangular */
/* >          matrix overwriting A (see below for further details). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D. */
/* >          If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >          interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* >          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and */
/* >          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* >          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) = */
/* >          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were */
/* >          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization */
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

/* > \ingroup doubleOTHERcomputational */

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
/* >  J. Lewis, Boeing Computer Services Company */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsptrf_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal t, r1, d11, d12, d21, d22;
    static integer kc, kk, kp;
    static doublereal wk;
    static integer kx, knc, kpc, npp;
    static doublereal wkm1, wkp1;
    static integer imax, jmax;
    extern /* Subroutine */ int dspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);
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

#line 206 "dsptrf.f"
    /* Parameter adjustments */
#line 206 "dsptrf.f"
    --ipiv;
#line 206 "dsptrf.f"
    --ap;
#line 206 "dsptrf.f"

#line 206 "dsptrf.f"
    /* Function Body */
#line 206 "dsptrf.f"
    *info = 0;
#line 207 "dsptrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 208 "dsptrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 209 "dsptrf.f"
	*info = -1;
#line 210 "dsptrf.f"
    } else if (*n < 0) {
#line 211 "dsptrf.f"
	*info = -2;
#line 212 "dsptrf.f"
    }
#line 213 "dsptrf.f"
    if (*info != 0) {
#line 214 "dsptrf.f"
	i__1 = -(*info);
#line 214 "dsptrf.f"
	xerbla_("DSPTRF", &i__1, (ftnlen)6);
#line 215 "dsptrf.f"
	return 0;
#line 216 "dsptrf.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 220 "dsptrf.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 222 "dsptrf.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 229 "dsptrf.f"
	k = *n;
#line 230 "dsptrf.f"
	kc = (*n - 1) * *n / 2 + 1;
#line 231 "dsptrf.f"
L10:
#line 232 "dsptrf.f"
	knc = kc;

/*        If K < 1, exit from loop */

#line 236 "dsptrf.f"
	if (k < 1) {
#line 236 "dsptrf.f"
	    goto L110;
#line 236 "dsptrf.f"
	}
#line 238 "dsptrf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 243 "dsptrf.f"
	absakk = (d__1 = ap[kc + k - 1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value */

#line 248 "dsptrf.f"
	if (k > 1) {
#line 249 "dsptrf.f"
	    i__1 = k - 1;
#line 249 "dsptrf.f"
	    imax = idamax_(&i__1, &ap[kc], &c__1);
#line 250 "dsptrf.f"
	    colmax = (d__1 = ap[kc + imax - 1], abs(d__1));
#line 251 "dsptrf.f"
	} else {
#line 252 "dsptrf.f"
	    colmax = 0.;
#line 253 "dsptrf.f"
	}

#line 255 "dsptrf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

#line 259 "dsptrf.f"
	    if (*info == 0) {
#line 259 "dsptrf.f"
		*info = k;
#line 259 "dsptrf.f"
	    }
#line 261 "dsptrf.f"
	    kp = k;
#line 262 "dsptrf.f"
	} else {
#line 263 "dsptrf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 267 "dsptrf.f"
		kp = k;
#line 268 "dsptrf.f"
	    } else {

#line 270 "dsptrf.f"
		rowmax = 0.;
#line 271 "dsptrf.f"
		jmax = imax;
#line 272 "dsptrf.f"
		kx = imax * (imax + 1) / 2 + imax;
#line 273 "dsptrf.f"
		i__1 = k;
#line 273 "dsptrf.f"
		for (j = imax + 1; j <= i__1; ++j) {
#line 274 "dsptrf.f"
		    if ((d__1 = ap[kx], abs(d__1)) > rowmax) {
#line 275 "dsptrf.f"
			rowmax = (d__1 = ap[kx], abs(d__1));
#line 276 "dsptrf.f"
			jmax = j;
#line 277 "dsptrf.f"
		    }
#line 278 "dsptrf.f"
		    kx += j;
#line 279 "dsptrf.f"
/* L20: */
#line 279 "dsptrf.f"
		}
#line 280 "dsptrf.f"
		kpc = (imax - 1) * imax / 2 + 1;
#line 281 "dsptrf.f"
		if (imax > 1) {
#line 282 "dsptrf.f"
		    i__1 = imax - 1;
#line 282 "dsptrf.f"
		    jmax = idamax_(&i__1, &ap[kpc], &c__1);
/* Computing MAX */
#line 283 "dsptrf.f"
		    d__2 = rowmax, d__3 = (d__1 = ap[kpc + jmax - 1], abs(
			    d__1));
#line 283 "dsptrf.f"
		    rowmax = max(d__2,d__3);
#line 284 "dsptrf.f"
		}

#line 286 "dsptrf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 290 "dsptrf.f"
		    kp = k;
#line 291 "dsptrf.f"
		} else if ((d__1 = ap[kpc + imax - 1], abs(d__1)) >= alpha * 
			rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 296 "dsptrf.f"
		    kp = imax;
#line 297 "dsptrf.f"
		} else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 302 "dsptrf.f"
		    kp = imax;
#line 303 "dsptrf.f"
		    kstep = 2;
#line 304 "dsptrf.f"
		}
#line 305 "dsptrf.f"
	    }

#line 307 "dsptrf.f"
	    kk = k - kstep + 1;
#line 308 "dsptrf.f"
	    if (kstep == 2) {
#line 308 "dsptrf.f"
		knc = knc - k + 1;
#line 308 "dsptrf.f"
	    }
#line 310 "dsptrf.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 315 "dsptrf.f"
		i__1 = kp - 1;
#line 315 "dsptrf.f"
		dswap_(&i__1, &ap[knc], &c__1, &ap[kpc], &c__1);
#line 316 "dsptrf.f"
		kx = kpc + kp - 1;
#line 317 "dsptrf.f"
		i__1 = kk - 1;
#line 317 "dsptrf.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 318 "dsptrf.f"
		    kx = kx + j - 1;
#line 319 "dsptrf.f"
		    t = ap[knc + j - 1];
#line 320 "dsptrf.f"
		    ap[knc + j - 1] = ap[kx];
#line 321 "dsptrf.f"
		    ap[kx] = t;
#line 322 "dsptrf.f"
/* L30: */
#line 322 "dsptrf.f"
		}
#line 323 "dsptrf.f"
		t = ap[knc + kk - 1];
#line 324 "dsptrf.f"
		ap[knc + kk - 1] = ap[kpc + kp - 1];
#line 325 "dsptrf.f"
		ap[kpc + kp - 1] = t;
#line 326 "dsptrf.f"
		if (kstep == 2) {
#line 327 "dsptrf.f"
		    t = ap[kc + k - 2];
#line 328 "dsptrf.f"
		    ap[kc + k - 2] = ap[kc + kp - 1];
#line 329 "dsptrf.f"
		    ap[kc + kp - 1] = t;
#line 330 "dsptrf.f"
		}
#line 331 "dsptrf.f"
	    }

/*           Update the leading submatrix */

#line 335 "dsptrf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */

#line 347 "dsptrf.f"
		r1 = 1. / ap[kc + k - 1];
#line 348 "dsptrf.f"
		i__1 = k - 1;
#line 348 "dsptrf.f"
		d__1 = -r1;
#line 348 "dsptrf.f"
		dspr_(uplo, &i__1, &d__1, &ap[kc], &c__1, &ap[1], (ftnlen)1);

/*              Store U(k) in column k */

#line 352 "dsptrf.f"
		i__1 = k - 1;
#line 352 "dsptrf.f"
		dscal_(&i__1, &r1, &ap[kc], &c__1);
#line 353 "dsptrf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */

#line 367 "dsptrf.f"
		if (k > 2) {

#line 369 "dsptrf.f"
		    d12 = ap[k - 1 + (k - 1) * k / 2];
#line 370 "dsptrf.f"
		    d22 = ap[k - 1 + (k - 2) * (k - 1) / 2] / d12;
#line 371 "dsptrf.f"
		    d11 = ap[k + (k - 1) * k / 2] / d12;
#line 372 "dsptrf.f"
		    t = 1. / (d11 * d22 - 1.);
#line 373 "dsptrf.f"
		    d12 = t / d12;

#line 375 "dsptrf.f"
		    for (j = k - 2; j >= 1; --j) {
#line 376 "dsptrf.f"
			wkm1 = d12 * (d11 * ap[j + (k - 2) * (k - 1) / 2] - 
				ap[j + (k - 1) * k / 2]);
#line 378 "dsptrf.f"
			wk = d12 * (d22 * ap[j + (k - 1) * k / 2] - ap[j + (k 
				- 2) * (k - 1) / 2]);
#line 380 "dsptrf.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 381 "dsptrf.f"
			    ap[i__ + (j - 1) * j / 2] = ap[i__ + (j - 1) * j /
				     2] - ap[i__ + (k - 1) * k / 2] * wk - ap[
				    i__ + (k - 2) * (k - 1) / 2] * wkm1;
#line 384 "dsptrf.f"
/* L40: */
#line 384 "dsptrf.f"
			}
#line 385 "dsptrf.f"
			ap[j + (k - 1) * k / 2] = wk;
#line 386 "dsptrf.f"
			ap[j + (k - 2) * (k - 1) / 2] = wkm1;
#line 387 "dsptrf.f"
/* L50: */
#line 387 "dsptrf.f"
		    }

#line 389 "dsptrf.f"
		}

#line 391 "dsptrf.f"
	    }
#line 392 "dsptrf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 396 "dsptrf.f"
	if (kstep == 1) {
#line 397 "dsptrf.f"
	    ipiv[k] = kp;
#line 398 "dsptrf.f"
	} else {
#line 399 "dsptrf.f"
	    ipiv[k] = -kp;
#line 400 "dsptrf.f"
	    ipiv[k - 1] = -kp;
#line 401 "dsptrf.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 405 "dsptrf.f"
	k -= kstep;
#line 406 "dsptrf.f"
	kc = knc - k;
#line 407 "dsptrf.f"
	goto L10;

#line 409 "dsptrf.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 416 "dsptrf.f"
	k = 1;
#line 417 "dsptrf.f"
	kc = 1;
#line 418 "dsptrf.f"
	npp = *n * (*n + 1) / 2;
#line 419 "dsptrf.f"
L60:
#line 420 "dsptrf.f"
	knc = kc;

/*        If K > N, exit from loop */

#line 424 "dsptrf.f"
	if (k > *n) {
#line 424 "dsptrf.f"
	    goto L110;
#line 424 "dsptrf.f"
	}
#line 426 "dsptrf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 431 "dsptrf.f"
	absakk = (d__1 = ap[kc], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value */

#line 436 "dsptrf.f"
	if (k < *n) {
#line 437 "dsptrf.f"
	    i__1 = *n - k;
#line 437 "dsptrf.f"
	    imax = k + idamax_(&i__1, &ap[kc + 1], &c__1);
#line 438 "dsptrf.f"
	    colmax = (d__1 = ap[kc + imax - k], abs(d__1));
#line 439 "dsptrf.f"
	} else {
#line 440 "dsptrf.f"
	    colmax = 0.;
#line 441 "dsptrf.f"
	}

#line 443 "dsptrf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

#line 447 "dsptrf.f"
	    if (*info == 0) {
#line 447 "dsptrf.f"
		*info = k;
#line 447 "dsptrf.f"
	    }
#line 449 "dsptrf.f"
	    kp = k;
#line 450 "dsptrf.f"
	} else {
#line 451 "dsptrf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 455 "dsptrf.f"
		kp = k;
#line 456 "dsptrf.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 461 "dsptrf.f"
		rowmax = 0.;
#line 462 "dsptrf.f"
		kx = kc + imax - k;
#line 463 "dsptrf.f"
		i__1 = imax - 1;
#line 463 "dsptrf.f"
		for (j = k; j <= i__1; ++j) {
#line 464 "dsptrf.f"
		    if ((d__1 = ap[kx], abs(d__1)) > rowmax) {
#line 465 "dsptrf.f"
			rowmax = (d__1 = ap[kx], abs(d__1));
#line 466 "dsptrf.f"
			jmax = j;
#line 467 "dsptrf.f"
		    }
#line 468 "dsptrf.f"
		    kx = kx + *n - j;
#line 469 "dsptrf.f"
/* L70: */
#line 469 "dsptrf.f"
		}
#line 470 "dsptrf.f"
		kpc = npp - (*n - imax + 1) * (*n - imax + 2) / 2 + 1;
#line 471 "dsptrf.f"
		if (imax < *n) {
#line 472 "dsptrf.f"
		    i__1 = *n - imax;
#line 472 "dsptrf.f"
		    jmax = imax + idamax_(&i__1, &ap[kpc + 1], &c__1);
/* Computing MAX */
#line 473 "dsptrf.f"
		    d__2 = rowmax, d__3 = (d__1 = ap[kpc + jmax - imax], abs(
			    d__1));
#line 473 "dsptrf.f"
		    rowmax = max(d__2,d__3);
#line 474 "dsptrf.f"
		}

#line 476 "dsptrf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 480 "dsptrf.f"
		    kp = k;
#line 481 "dsptrf.f"
		} else if ((d__1 = ap[kpc], abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 486 "dsptrf.f"
		    kp = imax;
#line 487 "dsptrf.f"
		} else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 492 "dsptrf.f"
		    kp = imax;
#line 493 "dsptrf.f"
		    kstep = 2;
#line 494 "dsptrf.f"
		}
#line 495 "dsptrf.f"
	    }

#line 497 "dsptrf.f"
	    kk = k + kstep - 1;
#line 498 "dsptrf.f"
	    if (kstep == 2) {
#line 498 "dsptrf.f"
		knc = knc + *n - k + 1;
#line 498 "dsptrf.f"
	    }
#line 500 "dsptrf.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 505 "dsptrf.f"
		if (kp < *n) {
#line 505 "dsptrf.f"
		    i__1 = *n - kp;
#line 505 "dsptrf.f"
		    dswap_(&i__1, &ap[knc + kp - kk + 1], &c__1, &ap[kpc + 1],
			     &c__1);
#line 505 "dsptrf.f"
		}
#line 508 "dsptrf.f"
		kx = knc + kp - kk;
#line 509 "dsptrf.f"
		i__1 = kp - 1;
#line 509 "dsptrf.f"
		for (j = kk + 1; j <= i__1; ++j) {
#line 510 "dsptrf.f"
		    kx = kx + *n - j + 1;
#line 511 "dsptrf.f"
		    t = ap[knc + j - kk];
#line 512 "dsptrf.f"
		    ap[knc + j - kk] = ap[kx];
#line 513 "dsptrf.f"
		    ap[kx] = t;
#line 514 "dsptrf.f"
/* L80: */
#line 514 "dsptrf.f"
		}
#line 515 "dsptrf.f"
		t = ap[knc];
#line 516 "dsptrf.f"
		ap[knc] = ap[kpc];
#line 517 "dsptrf.f"
		ap[kpc] = t;
#line 518 "dsptrf.f"
		if (kstep == 2) {
#line 519 "dsptrf.f"
		    t = ap[kc + 1];
#line 520 "dsptrf.f"
		    ap[kc + 1] = ap[kc + kp - k];
#line 521 "dsptrf.f"
		    ap[kc + kp - k] = t;
#line 522 "dsptrf.f"
		}
#line 523 "dsptrf.f"
	    }

/*           Update the trailing submatrix */

#line 527 "dsptrf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 535 "dsptrf.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */

#line 541 "dsptrf.f"
		    r1 = 1. / ap[kc];
#line 542 "dsptrf.f"
		    i__1 = *n - k;
#line 542 "dsptrf.f"
		    d__1 = -r1;
#line 542 "dsptrf.f"
		    dspr_(uplo, &i__1, &d__1, &ap[kc + 1], &c__1, &ap[kc + *n 
			    - k + 1], (ftnlen)1);

/*                 Store L(k) in column K */

#line 547 "dsptrf.f"
		    i__1 = *n - k;
#line 547 "dsptrf.f"
		    dscal_(&i__1, &r1, &ap[kc + 1], &c__1);
#line 548 "dsptrf.f"
		}
#line 549 "dsptrf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns K and K+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

#line 558 "dsptrf.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T */
/*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 568 "dsptrf.f"
		    d21 = ap[k + 1 + (k - 1) * ((*n << 1) - k) / 2];
#line 569 "dsptrf.f"
		    d11 = ap[k + 1 + k * ((*n << 1) - k - 1) / 2] / d21;
#line 570 "dsptrf.f"
		    d22 = ap[k + (k - 1) * ((*n << 1) - k) / 2] / d21;
#line 571 "dsptrf.f"
		    t = 1. / (d11 * d22 - 1.);
#line 572 "dsptrf.f"
		    d21 = t / d21;

#line 574 "dsptrf.f"
		    i__1 = *n;
#line 574 "dsptrf.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 575 "dsptrf.f"
			wk = d21 * (d11 * ap[j + (k - 1) * ((*n << 1) - k) / 
				2] - ap[j + k * ((*n << 1) - k - 1) / 2]);
#line 577 "dsptrf.f"
			wkp1 = d21 * (d22 * ap[j + k * ((*n << 1) - k - 1) / 
				2] - ap[j + (k - 1) * ((*n << 1) - k) / 2]);

#line 580 "dsptrf.f"
			i__2 = *n;
#line 580 "dsptrf.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 581 "dsptrf.f"
			    ap[i__ + (j - 1) * ((*n << 1) - j) / 2] = ap[i__ 
				    + (j - 1) * ((*n << 1) - j) / 2] - ap[i__ 
				    + (k - 1) * ((*n << 1) - k) / 2] * wk - 
				    ap[i__ + k * ((*n << 1) - k - 1) / 2] * 
				    wkp1;
#line 584 "dsptrf.f"
/* L90: */
#line 584 "dsptrf.f"
			}

#line 586 "dsptrf.f"
			ap[j + (k - 1) * ((*n << 1) - k) / 2] = wk;
#line 587 "dsptrf.f"
			ap[j + k * ((*n << 1) - k - 1) / 2] = wkp1;

#line 589 "dsptrf.f"
/* L100: */
#line 589 "dsptrf.f"
		    }
#line 590 "dsptrf.f"
		}
#line 591 "dsptrf.f"
	    }
#line 592 "dsptrf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 596 "dsptrf.f"
	if (kstep == 1) {
#line 597 "dsptrf.f"
	    ipiv[k] = kp;
#line 598 "dsptrf.f"
	} else {
#line 599 "dsptrf.f"
	    ipiv[k] = -kp;
#line 600 "dsptrf.f"
	    ipiv[k + 1] = -kp;
#line 601 "dsptrf.f"
	}

/*        Increase K and return to the start of the main loop */

#line 605 "dsptrf.f"
	k += kstep;
#line 606 "dsptrf.f"
	kc = knc + *n - k + 2;
#line 607 "dsptrf.f"
	goto L60;

#line 609 "dsptrf.f"
    }

#line 611 "dsptrf.f"
L110:
#line 612 "dsptrf.f"
    return 0;

/*     End of DSPTRF */

} /* dsptrf_ */

