#line 1 "ssptrf.f"
/* ssptrf.f -- translated by f2c (version 20100827).
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

#line 1 "ssptrf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SSPTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPTRF( UPLO, N, AP, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPTRF computes the factorization of a real symmetric matrix A stored */
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
/* >          AP is REAL array, dimension (N*(N+1)/2) */
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

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  5-96 - Based on modifications by J. Lewis, Boeing Computer Services */
/* >         Company */
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
/* > */
/*  ===================================================================== */
/* Subroutine */ int ssptrf_(char *uplo, integer *n, doublereal *ap, integer *
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
    extern /* Subroutine */ int sspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);
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
    static doublereal colmax, rowmax;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 204 "ssptrf.f"
    /* Parameter adjustments */
#line 204 "ssptrf.f"
    --ipiv;
#line 204 "ssptrf.f"
    --ap;
#line 204 "ssptrf.f"

#line 204 "ssptrf.f"
    /* Function Body */
#line 204 "ssptrf.f"
    *info = 0;
#line 205 "ssptrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 206 "ssptrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 207 "ssptrf.f"
	*info = -1;
#line 208 "ssptrf.f"
    } else if (*n < 0) {
#line 209 "ssptrf.f"
	*info = -2;
#line 210 "ssptrf.f"
    }
#line 211 "ssptrf.f"
    if (*info != 0) {
#line 212 "ssptrf.f"
	i__1 = -(*info);
#line 212 "ssptrf.f"
	xerbla_("SSPTRF", &i__1, (ftnlen)6);
#line 213 "ssptrf.f"
	return 0;
#line 214 "ssptrf.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 218 "ssptrf.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 220 "ssptrf.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 227 "ssptrf.f"
	k = *n;
#line 228 "ssptrf.f"
	kc = (*n - 1) * *n / 2 + 1;
#line 229 "ssptrf.f"
L10:
#line 230 "ssptrf.f"
	knc = kc;

/*        If K < 1, exit from loop */

#line 234 "ssptrf.f"
	if (k < 1) {
#line 234 "ssptrf.f"
	    goto L110;
#line 234 "ssptrf.f"
	}
#line 236 "ssptrf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 241 "ssptrf.f"
	absakk = (d__1 = ap[kc + k - 1], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value */

#line 246 "ssptrf.f"
	if (k > 1) {
#line 247 "ssptrf.f"
	    i__1 = k - 1;
#line 247 "ssptrf.f"
	    imax = isamax_(&i__1, &ap[kc], &c__1);
#line 248 "ssptrf.f"
	    colmax = (d__1 = ap[kc + imax - 1], abs(d__1));
#line 249 "ssptrf.f"
	} else {
#line 250 "ssptrf.f"
	    colmax = 0.;
#line 251 "ssptrf.f"
	}

#line 253 "ssptrf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

#line 257 "ssptrf.f"
	    if (*info == 0) {
#line 257 "ssptrf.f"
		*info = k;
#line 257 "ssptrf.f"
	    }
#line 259 "ssptrf.f"
	    kp = k;
#line 260 "ssptrf.f"
	} else {
#line 261 "ssptrf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 265 "ssptrf.f"
		kp = k;
#line 266 "ssptrf.f"
	    } else {

#line 268 "ssptrf.f"
		rowmax = 0.;
#line 269 "ssptrf.f"
		jmax = imax;
#line 270 "ssptrf.f"
		kx = imax * (imax + 1) / 2 + imax;
#line 271 "ssptrf.f"
		i__1 = k;
#line 271 "ssptrf.f"
		for (j = imax + 1; j <= i__1; ++j) {
#line 272 "ssptrf.f"
		    if ((d__1 = ap[kx], abs(d__1)) > rowmax) {
#line 273 "ssptrf.f"
			rowmax = (d__1 = ap[kx], abs(d__1));
#line 274 "ssptrf.f"
			jmax = j;
#line 275 "ssptrf.f"
		    }
#line 276 "ssptrf.f"
		    kx += j;
#line 277 "ssptrf.f"
/* L20: */
#line 277 "ssptrf.f"
		}
#line 278 "ssptrf.f"
		kpc = (imax - 1) * imax / 2 + 1;
#line 279 "ssptrf.f"
		if (imax > 1) {
#line 280 "ssptrf.f"
		    i__1 = imax - 1;
#line 280 "ssptrf.f"
		    jmax = isamax_(&i__1, &ap[kpc], &c__1);
/* Computing MAX */
#line 281 "ssptrf.f"
		    d__2 = rowmax, d__3 = (d__1 = ap[kpc + jmax - 1], abs(
			    d__1));
#line 281 "ssptrf.f"
		    rowmax = max(d__2,d__3);
#line 282 "ssptrf.f"
		}

#line 284 "ssptrf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 288 "ssptrf.f"
		    kp = k;
#line 289 "ssptrf.f"
		} else if ((d__1 = ap[kpc + imax - 1], abs(d__1)) >= alpha * 
			rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 294 "ssptrf.f"
		    kp = imax;
#line 295 "ssptrf.f"
		} else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 300 "ssptrf.f"
		    kp = imax;
#line 301 "ssptrf.f"
		    kstep = 2;
#line 302 "ssptrf.f"
		}
#line 303 "ssptrf.f"
	    }

#line 305 "ssptrf.f"
	    kk = k - kstep + 1;
#line 306 "ssptrf.f"
	    if (kstep == 2) {
#line 306 "ssptrf.f"
		knc = knc - k + 1;
#line 306 "ssptrf.f"
	    }
#line 308 "ssptrf.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 313 "ssptrf.f"
		i__1 = kp - 1;
#line 313 "ssptrf.f"
		sswap_(&i__1, &ap[knc], &c__1, &ap[kpc], &c__1);
#line 314 "ssptrf.f"
		kx = kpc + kp - 1;
#line 315 "ssptrf.f"
		i__1 = kk - 1;
#line 315 "ssptrf.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 316 "ssptrf.f"
		    kx = kx + j - 1;
#line 317 "ssptrf.f"
		    t = ap[knc + j - 1];
#line 318 "ssptrf.f"
		    ap[knc + j - 1] = ap[kx];
#line 319 "ssptrf.f"
		    ap[kx] = t;
#line 320 "ssptrf.f"
/* L30: */
#line 320 "ssptrf.f"
		}
#line 321 "ssptrf.f"
		t = ap[knc + kk - 1];
#line 322 "ssptrf.f"
		ap[knc + kk - 1] = ap[kpc + kp - 1];
#line 323 "ssptrf.f"
		ap[kpc + kp - 1] = t;
#line 324 "ssptrf.f"
		if (kstep == 2) {
#line 325 "ssptrf.f"
		    t = ap[kc + k - 2];
#line 326 "ssptrf.f"
		    ap[kc + k - 2] = ap[kc + kp - 1];
#line 327 "ssptrf.f"
		    ap[kc + kp - 1] = t;
#line 328 "ssptrf.f"
		}
#line 329 "ssptrf.f"
	    }

/*           Update the leading submatrix */

#line 333 "ssptrf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */

#line 345 "ssptrf.f"
		r1 = 1. / ap[kc + k - 1];
#line 346 "ssptrf.f"
		i__1 = k - 1;
#line 346 "ssptrf.f"
		d__1 = -r1;
#line 346 "ssptrf.f"
		sspr_(uplo, &i__1, &d__1, &ap[kc], &c__1, &ap[1], (ftnlen)1);

/*              Store U(k) in column k */

#line 350 "ssptrf.f"
		i__1 = k - 1;
#line 350 "ssptrf.f"
		sscal_(&i__1, &r1, &ap[kc], &c__1);
#line 351 "ssptrf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */

#line 365 "ssptrf.f"
		if (k > 2) {

#line 367 "ssptrf.f"
		    d12 = ap[k - 1 + (k - 1) * k / 2];
#line 368 "ssptrf.f"
		    d22 = ap[k - 1 + (k - 2) * (k - 1) / 2] / d12;
#line 369 "ssptrf.f"
		    d11 = ap[k + (k - 1) * k / 2] / d12;
#line 370 "ssptrf.f"
		    t = 1. / (d11 * d22 - 1.);
#line 371 "ssptrf.f"
		    d12 = t / d12;

#line 373 "ssptrf.f"
		    for (j = k - 2; j >= 1; --j) {
#line 374 "ssptrf.f"
			wkm1 = d12 * (d11 * ap[j + (k - 2) * (k - 1) / 2] - 
				ap[j + (k - 1) * k / 2]);
#line 376 "ssptrf.f"
			wk = d12 * (d22 * ap[j + (k - 1) * k / 2] - ap[j + (k 
				- 2) * (k - 1) / 2]);
#line 378 "ssptrf.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 379 "ssptrf.f"
			    ap[i__ + (j - 1) * j / 2] = ap[i__ + (j - 1) * j /
				     2] - ap[i__ + (k - 1) * k / 2] * wk - ap[
				    i__ + (k - 2) * (k - 1) / 2] * wkm1;
#line 382 "ssptrf.f"
/* L40: */
#line 382 "ssptrf.f"
			}
#line 383 "ssptrf.f"
			ap[j + (k - 1) * k / 2] = wk;
#line 384 "ssptrf.f"
			ap[j + (k - 2) * (k - 1) / 2] = wkm1;
#line 385 "ssptrf.f"
/* L50: */
#line 385 "ssptrf.f"
		    }

#line 387 "ssptrf.f"
		}

#line 389 "ssptrf.f"
	    }
#line 390 "ssptrf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 394 "ssptrf.f"
	if (kstep == 1) {
#line 395 "ssptrf.f"
	    ipiv[k] = kp;
#line 396 "ssptrf.f"
	} else {
#line 397 "ssptrf.f"
	    ipiv[k] = -kp;
#line 398 "ssptrf.f"
	    ipiv[k - 1] = -kp;
#line 399 "ssptrf.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 403 "ssptrf.f"
	k -= kstep;
#line 404 "ssptrf.f"
	kc = knc - k;
#line 405 "ssptrf.f"
	goto L10;

#line 407 "ssptrf.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 414 "ssptrf.f"
	k = 1;
#line 415 "ssptrf.f"
	kc = 1;
#line 416 "ssptrf.f"
	npp = *n * (*n + 1) / 2;
#line 417 "ssptrf.f"
L60:
#line 418 "ssptrf.f"
	knc = kc;

/*        If K > N, exit from loop */

#line 422 "ssptrf.f"
	if (k > *n) {
#line 422 "ssptrf.f"
	    goto L110;
#line 422 "ssptrf.f"
	}
#line 424 "ssptrf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 429 "ssptrf.f"
	absakk = (d__1 = ap[kc], abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value */

#line 434 "ssptrf.f"
	if (k < *n) {
#line 435 "ssptrf.f"
	    i__1 = *n - k;
#line 435 "ssptrf.f"
	    imax = k + isamax_(&i__1, &ap[kc + 1], &c__1);
#line 436 "ssptrf.f"
	    colmax = (d__1 = ap[kc + imax - k], abs(d__1));
#line 437 "ssptrf.f"
	} else {
#line 438 "ssptrf.f"
	    colmax = 0.;
#line 439 "ssptrf.f"
	}

#line 441 "ssptrf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

#line 445 "ssptrf.f"
	    if (*info == 0) {
#line 445 "ssptrf.f"
		*info = k;
#line 445 "ssptrf.f"
	    }
#line 447 "ssptrf.f"
	    kp = k;
#line 448 "ssptrf.f"
	} else {
#line 449 "ssptrf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 453 "ssptrf.f"
		kp = k;
#line 454 "ssptrf.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 459 "ssptrf.f"
		rowmax = 0.;
#line 460 "ssptrf.f"
		kx = kc + imax - k;
#line 461 "ssptrf.f"
		i__1 = imax - 1;
#line 461 "ssptrf.f"
		for (j = k; j <= i__1; ++j) {
#line 462 "ssptrf.f"
		    if ((d__1 = ap[kx], abs(d__1)) > rowmax) {
#line 463 "ssptrf.f"
			rowmax = (d__1 = ap[kx], abs(d__1));
#line 464 "ssptrf.f"
			jmax = j;
#line 465 "ssptrf.f"
		    }
#line 466 "ssptrf.f"
		    kx = kx + *n - j;
#line 467 "ssptrf.f"
/* L70: */
#line 467 "ssptrf.f"
		}
#line 468 "ssptrf.f"
		kpc = npp - (*n - imax + 1) * (*n - imax + 2) / 2 + 1;
#line 469 "ssptrf.f"
		if (imax < *n) {
#line 470 "ssptrf.f"
		    i__1 = *n - imax;
#line 470 "ssptrf.f"
		    jmax = imax + isamax_(&i__1, &ap[kpc + 1], &c__1);
/* Computing MAX */
#line 471 "ssptrf.f"
		    d__2 = rowmax, d__3 = (d__1 = ap[kpc + jmax - imax], abs(
			    d__1));
#line 471 "ssptrf.f"
		    rowmax = max(d__2,d__3);
#line 472 "ssptrf.f"
		}

#line 474 "ssptrf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 478 "ssptrf.f"
		    kp = k;
#line 479 "ssptrf.f"
		} else if ((d__1 = ap[kpc], abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 484 "ssptrf.f"
		    kp = imax;
#line 485 "ssptrf.f"
		} else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 490 "ssptrf.f"
		    kp = imax;
#line 491 "ssptrf.f"
		    kstep = 2;
#line 492 "ssptrf.f"
		}
#line 493 "ssptrf.f"
	    }

#line 495 "ssptrf.f"
	    kk = k + kstep - 1;
#line 496 "ssptrf.f"
	    if (kstep == 2) {
#line 496 "ssptrf.f"
		knc = knc + *n - k + 1;
#line 496 "ssptrf.f"
	    }
#line 498 "ssptrf.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 503 "ssptrf.f"
		if (kp < *n) {
#line 503 "ssptrf.f"
		    i__1 = *n - kp;
#line 503 "ssptrf.f"
		    sswap_(&i__1, &ap[knc + kp - kk + 1], &c__1, &ap[kpc + 1],
			     &c__1);
#line 503 "ssptrf.f"
		}
#line 506 "ssptrf.f"
		kx = knc + kp - kk;
#line 507 "ssptrf.f"
		i__1 = kp - 1;
#line 507 "ssptrf.f"
		for (j = kk + 1; j <= i__1; ++j) {
#line 508 "ssptrf.f"
		    kx = kx + *n - j + 1;
#line 509 "ssptrf.f"
		    t = ap[knc + j - kk];
#line 510 "ssptrf.f"
		    ap[knc + j - kk] = ap[kx];
#line 511 "ssptrf.f"
		    ap[kx] = t;
#line 512 "ssptrf.f"
/* L80: */
#line 512 "ssptrf.f"
		}
#line 513 "ssptrf.f"
		t = ap[knc];
#line 514 "ssptrf.f"
		ap[knc] = ap[kpc];
#line 515 "ssptrf.f"
		ap[kpc] = t;
#line 516 "ssptrf.f"
		if (kstep == 2) {
#line 517 "ssptrf.f"
		    t = ap[kc + 1];
#line 518 "ssptrf.f"
		    ap[kc + 1] = ap[kc + kp - k];
#line 519 "ssptrf.f"
		    ap[kc + kp - k] = t;
#line 520 "ssptrf.f"
		}
#line 521 "ssptrf.f"
	    }

/*           Update the trailing submatrix */

#line 525 "ssptrf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 533 "ssptrf.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */

#line 539 "ssptrf.f"
		    r1 = 1. / ap[kc];
#line 540 "ssptrf.f"
		    i__1 = *n - k;
#line 540 "ssptrf.f"
		    d__1 = -r1;
#line 540 "ssptrf.f"
		    sspr_(uplo, &i__1, &d__1, &ap[kc + 1], &c__1, &ap[kc + *n 
			    - k + 1], (ftnlen)1);

/*                 Store L(k) in column K */

#line 545 "ssptrf.f"
		    i__1 = *n - k;
#line 545 "ssptrf.f"
		    sscal_(&i__1, &r1, &ap[kc + 1], &c__1);
#line 546 "ssptrf.f"
		}
#line 547 "ssptrf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns K and K+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

#line 556 "ssptrf.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T */
/*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 566 "ssptrf.f"
		    d21 = ap[k + 1 + (k - 1) * ((*n << 1) - k) / 2];
#line 567 "ssptrf.f"
		    d11 = ap[k + 1 + k * ((*n << 1) - k - 1) / 2] / d21;
#line 568 "ssptrf.f"
		    d22 = ap[k + (k - 1) * ((*n << 1) - k) / 2] / d21;
#line 569 "ssptrf.f"
		    t = 1. / (d11 * d22 - 1.);
#line 570 "ssptrf.f"
		    d21 = t / d21;

#line 572 "ssptrf.f"
		    i__1 = *n;
#line 572 "ssptrf.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 573 "ssptrf.f"
			wk = d21 * (d11 * ap[j + (k - 1) * ((*n << 1) - k) / 
				2] - ap[j + k * ((*n << 1) - k - 1) / 2]);
#line 575 "ssptrf.f"
			wkp1 = d21 * (d22 * ap[j + k * ((*n << 1) - k - 1) / 
				2] - ap[j + (k - 1) * ((*n << 1) - k) / 2]);

#line 578 "ssptrf.f"
			i__2 = *n;
#line 578 "ssptrf.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 579 "ssptrf.f"
			    ap[i__ + (j - 1) * ((*n << 1) - j) / 2] = ap[i__ 
				    + (j - 1) * ((*n << 1) - j) / 2] - ap[i__ 
				    + (k - 1) * ((*n << 1) - k) / 2] * wk - 
				    ap[i__ + k * ((*n << 1) - k - 1) / 2] * 
				    wkp1;
#line 582 "ssptrf.f"
/* L90: */
#line 582 "ssptrf.f"
			}

#line 584 "ssptrf.f"
			ap[j + (k - 1) * ((*n << 1) - k) / 2] = wk;
#line 585 "ssptrf.f"
			ap[j + k * ((*n << 1) - k - 1) / 2] = wkp1;

#line 587 "ssptrf.f"
/* L100: */
#line 587 "ssptrf.f"
		    }
#line 588 "ssptrf.f"
		}
#line 589 "ssptrf.f"
	    }
#line 590 "ssptrf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 594 "ssptrf.f"
	if (kstep == 1) {
#line 595 "ssptrf.f"
	    ipiv[k] = kp;
#line 596 "ssptrf.f"
	} else {
#line 597 "ssptrf.f"
	    ipiv[k] = -kp;
#line 598 "ssptrf.f"
	    ipiv[k + 1] = -kp;
#line 599 "ssptrf.f"
	}

/*        Increase K and return to the start of the main loop */

#line 603 "ssptrf.f"
	k += kstep;
#line 604 "ssptrf.f"
	kc = knc + *n - k + 2;
#line 605 "ssptrf.f"
	goto L60;

#line 607 "ssptrf.f"
    }

#line 609 "ssptrf.f"
L110:
#line 610 "ssptrf.f"
    return 0;

/*     End of SSPTRF */

} /* ssptrf_ */

