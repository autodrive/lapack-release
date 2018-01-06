#line 1 "zhptrf.f"
/* zhptrf.f -- translated by f2c (version 20100827).
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

#line 1 "zhptrf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZHPTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPTRF( UPLO, N, AP, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPTRF computes the factorization of a complex Hermitian packed */
/* > matrix A using the Bunch-Kaufman diagonal pivoting method: */
/* > */
/* >    A = U*D*U**H  or  A = L*D*L**H */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is Hermitian and block diagonal with */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
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

/* > \ingroup complex16OTHERcomputational */

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
/* >  J. Lewis, Boeing Computer Services Company */

/*  ===================================================================== */
/* Subroutine */ int zhptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;
    static doublecomplex t;
    static doublereal r1, d11;
    static doublecomplex d12;
    static doublereal d22;
    static doublecomplex d21;
    static integer kc, kk, kp;
    static doublecomplex wk;
    static integer kx;
    static doublereal tt;
    static integer knc, kpc, npp;
    static doublecomplex wkm1, wkp1;
    static integer imax, jmax;
    extern /* Subroutine */ int zhpr_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal absakk;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal colmax;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static doublereal rowmax;


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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 214 "zhptrf.f"
    /* Parameter adjustments */
#line 214 "zhptrf.f"
    --ipiv;
#line 214 "zhptrf.f"
    --ap;
#line 214 "zhptrf.f"

#line 214 "zhptrf.f"
    /* Function Body */
#line 214 "zhptrf.f"
    *info = 0;
#line 215 "zhptrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 216 "zhptrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 217 "zhptrf.f"
	*info = -1;
#line 218 "zhptrf.f"
    } else if (*n < 0) {
#line 219 "zhptrf.f"
	*info = -2;
#line 220 "zhptrf.f"
    }
#line 221 "zhptrf.f"
    if (*info != 0) {
#line 222 "zhptrf.f"
	i__1 = -(*info);
#line 222 "zhptrf.f"
	xerbla_("ZHPTRF", &i__1, (ftnlen)6);
#line 223 "zhptrf.f"
	return 0;
#line 224 "zhptrf.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 228 "zhptrf.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 230 "zhptrf.f"
    if (upper) {

/*        Factorize A as U*D*U**H using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 237 "zhptrf.f"
	k = *n;
#line 238 "zhptrf.f"
	kc = (*n - 1) * *n / 2 + 1;
#line 239 "zhptrf.f"
L10:
#line 240 "zhptrf.f"
	knc = kc;

/*        If K < 1, exit from loop */

#line 244 "zhptrf.f"
	if (k < 1) {
#line 244 "zhptrf.f"
	    goto L110;
#line 244 "zhptrf.f"
	}
#line 246 "zhptrf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 251 "zhptrf.f"
	i__1 = kc + k - 1;
#line 251 "zhptrf.f"
	absakk = (d__1 = ap[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value */

#line 256 "zhptrf.f"
	if (k > 1) {
#line 257 "zhptrf.f"
	    i__1 = k - 1;
#line 257 "zhptrf.f"
	    imax = izamax_(&i__1, &ap[kc], &c__1);
#line 258 "zhptrf.f"
	    i__1 = kc + imax - 1;
#line 258 "zhptrf.f"
	    colmax = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kc + 
		    imax - 1]), abs(d__2));
#line 259 "zhptrf.f"
	} else {
#line 260 "zhptrf.f"
	    colmax = 0.;
#line 261 "zhptrf.f"
	}

#line 263 "zhptrf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

#line 267 "zhptrf.f"
	    if (*info == 0) {
#line 267 "zhptrf.f"
		*info = k;
#line 267 "zhptrf.f"
	    }
#line 269 "zhptrf.f"
	    kp = k;
#line 270 "zhptrf.f"
	    i__1 = kc + k - 1;
#line 270 "zhptrf.f"
	    i__2 = kc + k - 1;
#line 270 "zhptrf.f"
	    d__1 = ap[i__2].r;
#line 270 "zhptrf.f"
	    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 271 "zhptrf.f"
	} else {
#line 272 "zhptrf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 276 "zhptrf.f"
		kp = k;
#line 277 "zhptrf.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 282 "zhptrf.f"
		rowmax = 0.;
#line 283 "zhptrf.f"
		jmax = imax;
#line 284 "zhptrf.f"
		kx = imax * (imax + 1) / 2 + imax;
#line 285 "zhptrf.f"
		i__1 = k;
#line 285 "zhptrf.f"
		for (j = imax + 1; j <= i__1; ++j) {
#line 286 "zhptrf.f"
		    i__2 = kx;
#line 286 "zhptrf.f"
		    if ((d__1 = ap[i__2].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    kx]), abs(d__2)) > rowmax) {
#line 287 "zhptrf.f"
			i__2 = kx;
#line 287 "zhptrf.f"
			rowmax = (d__1 = ap[i__2].r, abs(d__1)) + (d__2 = 
				d_imag(&ap[kx]), abs(d__2));
#line 288 "zhptrf.f"
			jmax = j;
#line 289 "zhptrf.f"
		    }
#line 290 "zhptrf.f"
		    kx += j;
#line 291 "zhptrf.f"
/* L20: */
#line 291 "zhptrf.f"
		}
#line 292 "zhptrf.f"
		kpc = (imax - 1) * imax / 2 + 1;
#line 293 "zhptrf.f"
		if (imax > 1) {
#line 294 "zhptrf.f"
		    i__1 = imax - 1;
#line 294 "zhptrf.f"
		    jmax = izamax_(&i__1, &ap[kpc], &c__1);
/* Computing MAX */
#line 295 "zhptrf.f"
		    i__1 = kpc + jmax - 1;
#line 295 "zhptrf.f"
		    d__3 = rowmax, d__4 = (d__1 = ap[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&ap[kpc + jmax - 1]), abs(d__2));
#line 295 "zhptrf.f"
		    rowmax = max(d__3,d__4);
#line 296 "zhptrf.f"
		}

#line 298 "zhptrf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 302 "zhptrf.f"
		    kp = k;
#line 303 "zhptrf.f"
		} else /* if(complicated condition) */ {
#line 303 "zhptrf.f"
		    i__1 = kpc + imax - 1;
#line 303 "zhptrf.f"
		    if ((d__1 = ap[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 309 "zhptrf.f"
			kp = imax;
#line 310 "zhptrf.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 315 "zhptrf.f"
			kp = imax;
#line 316 "zhptrf.f"
			kstep = 2;
#line 317 "zhptrf.f"
		    }
#line 317 "zhptrf.f"
		}
#line 318 "zhptrf.f"
	    }

#line 320 "zhptrf.f"
	    kk = k - kstep + 1;
#line 321 "zhptrf.f"
	    if (kstep == 2) {
#line 321 "zhptrf.f"
		knc = knc - k + 1;
#line 321 "zhptrf.f"
	    }
#line 323 "zhptrf.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 328 "zhptrf.f"
		i__1 = kp - 1;
#line 328 "zhptrf.f"
		zswap_(&i__1, &ap[knc], &c__1, &ap[kpc], &c__1);
#line 329 "zhptrf.f"
		kx = kpc + kp - 1;
#line 330 "zhptrf.f"
		i__1 = kk - 1;
#line 330 "zhptrf.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 331 "zhptrf.f"
		    kx = kx + j - 1;
#line 332 "zhptrf.f"
		    d_cnjg(&z__1, &ap[knc + j - 1]);
#line 332 "zhptrf.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 333 "zhptrf.f"
		    i__2 = knc + j - 1;
#line 333 "zhptrf.f"
		    d_cnjg(&z__1, &ap[kx]);
#line 333 "zhptrf.f"
		    ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 334 "zhptrf.f"
		    i__2 = kx;
#line 334 "zhptrf.f"
		    ap[i__2].r = t.r, ap[i__2].i = t.i;
#line 335 "zhptrf.f"
/* L30: */
#line 335 "zhptrf.f"
		}
#line 336 "zhptrf.f"
		i__1 = kx + kk - 1;
#line 336 "zhptrf.f"
		d_cnjg(&z__1, &ap[kx + kk - 1]);
#line 336 "zhptrf.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 337 "zhptrf.f"
		i__1 = knc + kk - 1;
#line 337 "zhptrf.f"
		r1 = ap[i__1].r;
#line 338 "zhptrf.f"
		i__1 = knc + kk - 1;
#line 338 "zhptrf.f"
		i__2 = kpc + kp - 1;
#line 338 "zhptrf.f"
		d__1 = ap[i__2].r;
#line 338 "zhptrf.f"
		ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 339 "zhptrf.f"
		i__1 = kpc + kp - 1;
#line 339 "zhptrf.f"
		ap[i__1].r = r1, ap[i__1].i = 0.;
#line 340 "zhptrf.f"
		if (kstep == 2) {
#line 341 "zhptrf.f"
		    i__1 = kc + k - 1;
#line 341 "zhptrf.f"
		    i__2 = kc + k - 1;
#line 341 "zhptrf.f"
		    d__1 = ap[i__2].r;
#line 341 "zhptrf.f"
		    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 342 "zhptrf.f"
		    i__1 = kc + k - 2;
#line 342 "zhptrf.f"
		    t.r = ap[i__1].r, t.i = ap[i__1].i;
#line 343 "zhptrf.f"
		    i__1 = kc + k - 2;
#line 343 "zhptrf.f"
		    i__2 = kc + kp - 1;
#line 343 "zhptrf.f"
		    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 344 "zhptrf.f"
		    i__1 = kc + kp - 1;
#line 344 "zhptrf.f"
		    ap[i__1].r = t.r, ap[i__1].i = t.i;
#line 345 "zhptrf.f"
		}
#line 346 "zhptrf.f"
	    } else {
#line 347 "zhptrf.f"
		i__1 = kc + k - 1;
#line 347 "zhptrf.f"
		i__2 = kc + k - 1;
#line 347 "zhptrf.f"
		d__1 = ap[i__2].r;
#line 347 "zhptrf.f"
		ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 348 "zhptrf.f"
		if (kstep == 2) {
#line 348 "zhptrf.f"
		    i__1 = kc - 1;
#line 348 "zhptrf.f"
		    i__2 = kc - 1;
#line 348 "zhptrf.f"
		    d__1 = ap[i__2].r;
#line 348 "zhptrf.f"
		    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 348 "zhptrf.f"
		}
#line 350 "zhptrf.f"
	    }

/*           Update the leading submatrix */

#line 354 "zhptrf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H */

#line 366 "zhptrf.f"
		i__1 = kc + k - 1;
#line 366 "zhptrf.f"
		r1 = 1. / ap[i__1].r;
#line 367 "zhptrf.f"
		i__1 = k - 1;
#line 367 "zhptrf.f"
		d__1 = -r1;
#line 367 "zhptrf.f"
		zhpr_(uplo, &i__1, &d__1, &ap[kc], &c__1, &ap[1], (ftnlen)1);

/*              Store U(k) in column k */

#line 371 "zhptrf.f"
		i__1 = k - 1;
#line 371 "zhptrf.f"
		zdscal_(&i__1, &r1, &ap[kc], &c__1);
#line 372 "zhptrf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H */

#line 386 "zhptrf.f"
		if (k > 2) {

#line 388 "zhptrf.f"
		    i__1 = k - 1 + (k - 1) * k / 2;
#line 388 "zhptrf.f"
		    d__1 = ap[i__1].r;
#line 388 "zhptrf.f"
		    d__2 = d_imag(&ap[k - 1 + (k - 1) * k / 2]);
#line 388 "zhptrf.f"
		    d__ = dlapy2_(&d__1, &d__2);
#line 390 "zhptrf.f"
		    i__1 = k - 1 + (k - 2) * (k - 1) / 2;
#line 390 "zhptrf.f"
		    d22 = ap[i__1].r / d__;
#line 391 "zhptrf.f"
		    i__1 = k + (k - 1) * k / 2;
#line 391 "zhptrf.f"
		    d11 = ap[i__1].r / d__;
#line 392 "zhptrf.f"
		    tt = 1. / (d11 * d22 - 1.);
#line 393 "zhptrf.f"
		    i__1 = k - 1 + (k - 1) * k / 2;
#line 393 "zhptrf.f"
		    z__1.r = ap[i__1].r / d__, z__1.i = ap[i__1].i / d__;
#line 393 "zhptrf.f"
		    d12.r = z__1.r, d12.i = z__1.i;
#line 394 "zhptrf.f"
		    d__ = tt / d__;

#line 396 "zhptrf.f"
		    for (j = k - 2; j >= 1; --j) {
#line 397 "zhptrf.f"
			i__1 = j + (k - 2) * (k - 1) / 2;
#line 397 "zhptrf.f"
			z__3.r = d11 * ap[i__1].r, z__3.i = d11 * ap[i__1].i;
#line 397 "zhptrf.f"
			d_cnjg(&z__5, &d12);
#line 397 "zhptrf.f"
			i__2 = j + (k - 1) * k / 2;
#line 397 "zhptrf.f"
			z__4.r = z__5.r * ap[i__2].r - z__5.i * ap[i__2].i, 
				z__4.i = z__5.r * ap[i__2].i + z__5.i * ap[
				i__2].r;
#line 397 "zhptrf.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 397 "zhptrf.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 397 "zhptrf.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 399 "zhptrf.f"
			i__1 = j + (k - 1) * k / 2;
#line 399 "zhptrf.f"
			z__3.r = d22 * ap[i__1].r, z__3.i = d22 * ap[i__1].i;
#line 399 "zhptrf.f"
			i__2 = j + (k - 2) * (k - 1) / 2;
#line 399 "zhptrf.f"
			z__4.r = d12.r * ap[i__2].r - d12.i * ap[i__2].i, 
				z__4.i = d12.r * ap[i__2].i + d12.i * ap[i__2]
				.r;
#line 399 "zhptrf.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 399 "zhptrf.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 399 "zhptrf.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 401 "zhptrf.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 402 "zhptrf.f"
			    i__1 = i__ + (j - 1) * j / 2;
#line 402 "zhptrf.f"
			    i__2 = i__ + (j - 1) * j / 2;
#line 402 "zhptrf.f"
			    i__3 = i__ + (k - 1) * k / 2;
#line 402 "zhptrf.f"
			    d_cnjg(&z__4, &wk);
#line 402 "zhptrf.f"
			    z__3.r = ap[i__3].r * z__4.r - ap[i__3].i * 
				    z__4.i, z__3.i = ap[i__3].r * z__4.i + ap[
				    i__3].i * z__4.r;
#line 402 "zhptrf.f"
			    z__2.r = ap[i__2].r - z__3.r, z__2.i = ap[i__2].i 
				    - z__3.i;
#line 402 "zhptrf.f"
			    i__4 = i__ + (k - 2) * (k - 1) / 2;
#line 402 "zhptrf.f"
			    d_cnjg(&z__6, &wkm1);
#line 402 "zhptrf.f"
			    z__5.r = ap[i__4].r * z__6.r - ap[i__4].i * 
				    z__6.i, z__5.i = ap[i__4].r * z__6.i + ap[
				    i__4].i * z__6.r;
#line 402 "zhptrf.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 402 "zhptrf.f"
			    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 405 "zhptrf.f"
/* L40: */
#line 405 "zhptrf.f"
			}
#line 406 "zhptrf.f"
			i__1 = j + (k - 1) * k / 2;
#line 406 "zhptrf.f"
			ap[i__1].r = wk.r, ap[i__1].i = wk.i;
#line 407 "zhptrf.f"
			i__1 = j + (k - 2) * (k - 1) / 2;
#line 407 "zhptrf.f"
			ap[i__1].r = wkm1.r, ap[i__1].i = wkm1.i;
#line 408 "zhptrf.f"
			i__1 = j + (j - 1) * j / 2;
#line 408 "zhptrf.f"
			i__2 = j + (j - 1) * j / 2;
#line 408 "zhptrf.f"
			d__1 = ap[i__2].r;
#line 408 "zhptrf.f"
			z__1.r = d__1, z__1.i = 0.;
#line 408 "zhptrf.f"
			ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 410 "zhptrf.f"
/* L50: */
#line 410 "zhptrf.f"
		    }

#line 412 "zhptrf.f"
		}

#line 414 "zhptrf.f"
	    }
#line 415 "zhptrf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 419 "zhptrf.f"
	if (kstep == 1) {
#line 420 "zhptrf.f"
	    ipiv[k] = kp;
#line 421 "zhptrf.f"
	} else {
#line 422 "zhptrf.f"
	    ipiv[k] = -kp;
#line 423 "zhptrf.f"
	    ipiv[k - 1] = -kp;
#line 424 "zhptrf.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 428 "zhptrf.f"
	k -= kstep;
#line 429 "zhptrf.f"
	kc = knc - k;
#line 430 "zhptrf.f"
	goto L10;

#line 432 "zhptrf.f"
    } else {

/*        Factorize A as L*D*L**H using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 439 "zhptrf.f"
	k = 1;
#line 440 "zhptrf.f"
	kc = 1;
#line 441 "zhptrf.f"
	npp = *n * (*n + 1) / 2;
#line 442 "zhptrf.f"
L60:
#line 443 "zhptrf.f"
	knc = kc;

/*        If K > N, exit from loop */

#line 447 "zhptrf.f"
	if (k > *n) {
#line 447 "zhptrf.f"
	    goto L110;
#line 447 "zhptrf.f"
	}
#line 449 "zhptrf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 454 "zhptrf.f"
	i__1 = kc;
#line 454 "zhptrf.f"
	absakk = (d__1 = ap[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value */

#line 459 "zhptrf.f"
	if (k < *n) {
#line 460 "zhptrf.f"
	    i__1 = *n - k;
#line 460 "zhptrf.f"
	    imax = k + izamax_(&i__1, &ap[kc + 1], &c__1);
#line 461 "zhptrf.f"
	    i__1 = kc + imax - k;
#line 461 "zhptrf.f"
	    colmax = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kc + 
		    imax - k]), abs(d__2));
#line 462 "zhptrf.f"
	} else {
#line 463 "zhptrf.f"
	    colmax = 0.;
#line 464 "zhptrf.f"
	}

#line 466 "zhptrf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

#line 470 "zhptrf.f"
	    if (*info == 0) {
#line 470 "zhptrf.f"
		*info = k;
#line 470 "zhptrf.f"
	    }
#line 472 "zhptrf.f"
	    kp = k;
#line 473 "zhptrf.f"
	    i__1 = kc;
#line 473 "zhptrf.f"
	    i__2 = kc;
#line 473 "zhptrf.f"
	    d__1 = ap[i__2].r;
#line 473 "zhptrf.f"
	    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 474 "zhptrf.f"
	} else {
#line 475 "zhptrf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 479 "zhptrf.f"
		kp = k;
#line 480 "zhptrf.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 485 "zhptrf.f"
		rowmax = 0.;
#line 486 "zhptrf.f"
		kx = kc + imax - k;
#line 487 "zhptrf.f"
		i__1 = imax - 1;
#line 487 "zhptrf.f"
		for (j = k; j <= i__1; ++j) {
#line 488 "zhptrf.f"
		    i__2 = kx;
#line 488 "zhptrf.f"
		    if ((d__1 = ap[i__2].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    kx]), abs(d__2)) > rowmax) {
#line 489 "zhptrf.f"
			i__2 = kx;
#line 489 "zhptrf.f"
			rowmax = (d__1 = ap[i__2].r, abs(d__1)) + (d__2 = 
				d_imag(&ap[kx]), abs(d__2));
#line 490 "zhptrf.f"
			jmax = j;
#line 491 "zhptrf.f"
		    }
#line 492 "zhptrf.f"
		    kx = kx + *n - j;
#line 493 "zhptrf.f"
/* L70: */
#line 493 "zhptrf.f"
		}
#line 494 "zhptrf.f"
		kpc = npp - (*n - imax + 1) * (*n - imax + 2) / 2 + 1;
#line 495 "zhptrf.f"
		if (imax < *n) {
#line 496 "zhptrf.f"
		    i__1 = *n - imax;
#line 496 "zhptrf.f"
		    jmax = imax + izamax_(&i__1, &ap[kpc + 1], &c__1);
/* Computing MAX */
#line 497 "zhptrf.f"
		    i__1 = kpc + jmax - imax;
#line 497 "zhptrf.f"
		    d__3 = rowmax, d__4 = (d__1 = ap[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&ap[kpc + jmax - imax]), abs(d__2));
#line 497 "zhptrf.f"
		    rowmax = max(d__3,d__4);
#line 498 "zhptrf.f"
		}

#line 500 "zhptrf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 504 "zhptrf.f"
		    kp = k;
#line 505 "zhptrf.f"
		} else /* if(complicated condition) */ {
#line 505 "zhptrf.f"
		    i__1 = kpc;
#line 505 "zhptrf.f"
		    if ((d__1 = ap[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 510 "zhptrf.f"
			kp = imax;
#line 511 "zhptrf.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 516 "zhptrf.f"
			kp = imax;
#line 517 "zhptrf.f"
			kstep = 2;
#line 518 "zhptrf.f"
		    }
#line 518 "zhptrf.f"
		}
#line 519 "zhptrf.f"
	    }

#line 521 "zhptrf.f"
	    kk = k + kstep - 1;
#line 522 "zhptrf.f"
	    if (kstep == 2) {
#line 522 "zhptrf.f"
		knc = knc + *n - k + 1;
#line 522 "zhptrf.f"
	    }
#line 524 "zhptrf.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 529 "zhptrf.f"
		if (kp < *n) {
#line 529 "zhptrf.f"
		    i__1 = *n - kp;
#line 529 "zhptrf.f"
		    zswap_(&i__1, &ap[knc + kp - kk + 1], &c__1, &ap[kpc + 1],
			     &c__1);
#line 529 "zhptrf.f"
		}
#line 532 "zhptrf.f"
		kx = knc + kp - kk;
#line 533 "zhptrf.f"
		i__1 = kp - 1;
#line 533 "zhptrf.f"
		for (j = kk + 1; j <= i__1; ++j) {
#line 534 "zhptrf.f"
		    kx = kx + *n - j + 1;
#line 535 "zhptrf.f"
		    d_cnjg(&z__1, &ap[knc + j - kk]);
#line 535 "zhptrf.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 536 "zhptrf.f"
		    i__2 = knc + j - kk;
#line 536 "zhptrf.f"
		    d_cnjg(&z__1, &ap[kx]);
#line 536 "zhptrf.f"
		    ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 537 "zhptrf.f"
		    i__2 = kx;
#line 537 "zhptrf.f"
		    ap[i__2].r = t.r, ap[i__2].i = t.i;
#line 538 "zhptrf.f"
/* L80: */
#line 538 "zhptrf.f"
		}
#line 539 "zhptrf.f"
		i__1 = knc + kp - kk;
#line 539 "zhptrf.f"
		d_cnjg(&z__1, &ap[knc + kp - kk]);
#line 539 "zhptrf.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 540 "zhptrf.f"
		i__1 = knc;
#line 540 "zhptrf.f"
		r1 = ap[i__1].r;
#line 541 "zhptrf.f"
		i__1 = knc;
#line 541 "zhptrf.f"
		i__2 = kpc;
#line 541 "zhptrf.f"
		d__1 = ap[i__2].r;
#line 541 "zhptrf.f"
		ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 542 "zhptrf.f"
		i__1 = kpc;
#line 542 "zhptrf.f"
		ap[i__1].r = r1, ap[i__1].i = 0.;
#line 543 "zhptrf.f"
		if (kstep == 2) {
#line 544 "zhptrf.f"
		    i__1 = kc;
#line 544 "zhptrf.f"
		    i__2 = kc;
#line 544 "zhptrf.f"
		    d__1 = ap[i__2].r;
#line 544 "zhptrf.f"
		    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 545 "zhptrf.f"
		    i__1 = kc + 1;
#line 545 "zhptrf.f"
		    t.r = ap[i__1].r, t.i = ap[i__1].i;
#line 546 "zhptrf.f"
		    i__1 = kc + 1;
#line 546 "zhptrf.f"
		    i__2 = kc + kp - k;
#line 546 "zhptrf.f"
		    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 547 "zhptrf.f"
		    i__1 = kc + kp - k;
#line 547 "zhptrf.f"
		    ap[i__1].r = t.r, ap[i__1].i = t.i;
#line 548 "zhptrf.f"
		}
#line 549 "zhptrf.f"
	    } else {
#line 550 "zhptrf.f"
		i__1 = kc;
#line 550 "zhptrf.f"
		i__2 = kc;
#line 550 "zhptrf.f"
		d__1 = ap[i__2].r;
#line 550 "zhptrf.f"
		ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 551 "zhptrf.f"
		if (kstep == 2) {
#line 551 "zhptrf.f"
		    i__1 = knc;
#line 551 "zhptrf.f"
		    i__2 = knc;
#line 551 "zhptrf.f"
		    d__1 = ap[i__2].r;
#line 551 "zhptrf.f"
		    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 551 "zhptrf.f"
		}
#line 553 "zhptrf.f"
	    }

/*           Update the trailing submatrix */

#line 557 "zhptrf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 565 "zhptrf.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H */

#line 571 "zhptrf.f"
		    i__1 = kc;
#line 571 "zhptrf.f"
		    r1 = 1. / ap[i__1].r;
#line 572 "zhptrf.f"
		    i__1 = *n - k;
#line 572 "zhptrf.f"
		    d__1 = -r1;
#line 572 "zhptrf.f"
		    zhpr_(uplo, &i__1, &d__1, &ap[kc + 1], &c__1, &ap[kc + *n 
			    - k + 1], (ftnlen)1);

/*                 Store L(k) in column K */

#line 577 "zhptrf.f"
		    i__1 = *n - k;
#line 577 "zhptrf.f"
		    zdscal_(&i__1, &r1, &ap[kc + 1], &c__1);
#line 578 "zhptrf.f"
		}
#line 579 "zhptrf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns K and K+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

#line 588 "zhptrf.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H */
/*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 598 "zhptrf.f"
		    i__1 = k + 1 + (k - 1) * ((*n << 1) - k) / 2;
#line 598 "zhptrf.f"
		    d__1 = ap[i__1].r;
#line 598 "zhptrf.f"
		    d__2 = d_imag(&ap[k + 1 + (k - 1) * ((*n << 1) - k) / 2]);
#line 598 "zhptrf.f"
		    d__ = dlapy2_(&d__1, &d__2);
#line 600 "zhptrf.f"
		    i__1 = k + 1 + k * ((*n << 1) - k - 1) / 2;
#line 600 "zhptrf.f"
		    d11 = ap[i__1].r / d__;
#line 601 "zhptrf.f"
		    i__1 = k + (k - 1) * ((*n << 1) - k) / 2;
#line 601 "zhptrf.f"
		    d22 = ap[i__1].r / d__;
#line 602 "zhptrf.f"
		    tt = 1. / (d11 * d22 - 1.);
#line 603 "zhptrf.f"
		    i__1 = k + 1 + (k - 1) * ((*n << 1) - k) / 2;
#line 603 "zhptrf.f"
		    z__1.r = ap[i__1].r / d__, z__1.i = ap[i__1].i / d__;
#line 603 "zhptrf.f"
		    d21.r = z__1.r, d21.i = z__1.i;
#line 604 "zhptrf.f"
		    d__ = tt / d__;

#line 606 "zhptrf.f"
		    i__1 = *n;
#line 606 "zhptrf.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 607 "zhptrf.f"
			i__2 = j + (k - 1) * ((*n << 1) - k) / 2;
#line 607 "zhptrf.f"
			z__3.r = d11 * ap[i__2].r, z__3.i = d11 * ap[i__2].i;
#line 607 "zhptrf.f"
			i__3 = j + k * ((*n << 1) - k - 1) / 2;
#line 607 "zhptrf.f"
			z__4.r = d21.r * ap[i__3].r - d21.i * ap[i__3].i, 
				z__4.i = d21.r * ap[i__3].i + d21.i * ap[i__3]
				.r;
#line 607 "zhptrf.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 607 "zhptrf.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 607 "zhptrf.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 609 "zhptrf.f"
			i__2 = j + k * ((*n << 1) - k - 1) / 2;
#line 609 "zhptrf.f"
			z__3.r = d22 * ap[i__2].r, z__3.i = d22 * ap[i__2].i;
#line 609 "zhptrf.f"
			d_cnjg(&z__5, &d21);
#line 609 "zhptrf.f"
			i__3 = j + (k - 1) * ((*n << 1) - k) / 2;
#line 609 "zhptrf.f"
			z__4.r = z__5.r * ap[i__3].r - z__5.i * ap[i__3].i, 
				z__4.i = z__5.r * ap[i__3].i + z__5.i * ap[
				i__3].r;
#line 609 "zhptrf.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 609 "zhptrf.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 609 "zhptrf.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;
#line 612 "zhptrf.f"
			i__2 = *n;
#line 612 "zhptrf.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 613 "zhptrf.f"
			    i__3 = i__ + (j - 1) * ((*n << 1) - j) / 2;
#line 613 "zhptrf.f"
			    i__4 = i__ + (j - 1) * ((*n << 1) - j) / 2;
#line 613 "zhptrf.f"
			    i__5 = i__ + (k - 1) * ((*n << 1) - k) / 2;
#line 613 "zhptrf.f"
			    d_cnjg(&z__4, &wk);
#line 613 "zhptrf.f"
			    z__3.r = ap[i__5].r * z__4.r - ap[i__5].i * 
				    z__4.i, z__3.i = ap[i__5].r * z__4.i + ap[
				    i__5].i * z__4.r;
#line 613 "zhptrf.f"
			    z__2.r = ap[i__4].r - z__3.r, z__2.i = ap[i__4].i 
				    - z__3.i;
#line 613 "zhptrf.f"
			    i__6 = i__ + k * ((*n << 1) - k - 1) / 2;
#line 613 "zhptrf.f"
			    d_cnjg(&z__6, &wkp1);
#line 613 "zhptrf.f"
			    z__5.r = ap[i__6].r * z__6.r - ap[i__6].i * 
				    z__6.i, z__5.i = ap[i__6].r * z__6.i + ap[
				    i__6].i * z__6.r;
#line 613 "zhptrf.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 613 "zhptrf.f"
			    ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 617 "zhptrf.f"
/* L90: */
#line 617 "zhptrf.f"
			}
#line 618 "zhptrf.f"
			i__2 = j + (k - 1) * ((*n << 1) - k) / 2;
#line 618 "zhptrf.f"
			ap[i__2].r = wk.r, ap[i__2].i = wk.i;
#line 619 "zhptrf.f"
			i__2 = j + k * ((*n << 1) - k - 1) / 2;
#line 619 "zhptrf.f"
			ap[i__2].r = wkp1.r, ap[i__2].i = wkp1.i;
#line 620 "zhptrf.f"
			i__2 = j + (j - 1) * ((*n << 1) - j) / 2;
#line 620 "zhptrf.f"
			i__3 = j + (j - 1) * ((*n << 1) - j) / 2;
#line 620 "zhptrf.f"
			d__1 = ap[i__3].r;
#line 620 "zhptrf.f"
			z__1.r = d__1, z__1.i = 0.;
#line 620 "zhptrf.f"
			ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 623 "zhptrf.f"
/* L100: */
#line 623 "zhptrf.f"
		    }
#line 624 "zhptrf.f"
		}
#line 625 "zhptrf.f"
	    }
#line 626 "zhptrf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 630 "zhptrf.f"
	if (kstep == 1) {
#line 631 "zhptrf.f"
	    ipiv[k] = kp;
#line 632 "zhptrf.f"
	} else {
#line 633 "zhptrf.f"
	    ipiv[k] = -kp;
#line 634 "zhptrf.f"
	    ipiv[k + 1] = -kp;
#line 635 "zhptrf.f"
	}

/*        Increase K and return to the start of the main loop */

#line 639 "zhptrf.f"
	k += kstep;
#line 640 "zhptrf.f"
	kc = knc + *n - k + 2;
#line 641 "zhptrf.f"
	goto L60;

#line 643 "zhptrf.f"
    }

#line 645 "zhptrf.f"
L110:
#line 646 "zhptrf.f"
    return 0;

/*     End of ZHPTRF */

} /* zhptrf_ */

