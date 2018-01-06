#line 1 "csptrf.f"
/* csptrf.f -- translated by f2c (version 20100827).
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

#line 1 "csptrf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CSPTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csptrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csptrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csptrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSPTRF( UPLO, N, AP, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSPTRF computes the factorization of a complex symmetric matrix A */
/* > stored in packed format using the Bunch-Kaufman diagonal pivoting */
/* > method: */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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

/* > \ingroup complexOTHERcomputational */

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
/* Subroutine */ int csptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex t, r1, d11, d12, d21, d22;
    static integer kc, kk, kp;
    static doublecomplex wk;
    static integer kx, knc, kpc, npp;
    static doublecomplex wkm1, wkp1;
    static integer imax, jmax;
    extern /* Subroutine */ int cspr_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer kstep;
    static logical upper;
    static doublereal absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 213 "csptrf.f"
    /* Parameter adjustments */
#line 213 "csptrf.f"
    --ipiv;
#line 213 "csptrf.f"
    --ap;
#line 213 "csptrf.f"

#line 213 "csptrf.f"
    /* Function Body */
#line 213 "csptrf.f"
    *info = 0;
#line 214 "csptrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 215 "csptrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 216 "csptrf.f"
	*info = -1;
#line 217 "csptrf.f"
    } else if (*n < 0) {
#line 218 "csptrf.f"
	*info = -2;
#line 219 "csptrf.f"
    }
#line 220 "csptrf.f"
    if (*info != 0) {
#line 221 "csptrf.f"
	i__1 = -(*info);
#line 221 "csptrf.f"
	xerbla_("CSPTRF", &i__1, (ftnlen)6);
#line 222 "csptrf.f"
	return 0;
#line 223 "csptrf.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 227 "csptrf.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 229 "csptrf.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 236 "csptrf.f"
	k = *n;
#line 237 "csptrf.f"
	kc = (*n - 1) * *n / 2 + 1;
#line 238 "csptrf.f"
L10:
#line 239 "csptrf.f"
	knc = kc;

/*        If K < 1, exit from loop */

#line 243 "csptrf.f"
	if (k < 1) {
#line 243 "csptrf.f"
	    goto L110;
#line 243 "csptrf.f"
	}
#line 245 "csptrf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 250 "csptrf.f"
	i__1 = kc + k - 1;
#line 250 "csptrf.f"
	absakk = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kc + k - 
		1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value */

#line 255 "csptrf.f"
	if (k > 1) {
#line 256 "csptrf.f"
	    i__1 = k - 1;
#line 256 "csptrf.f"
	    imax = icamax_(&i__1, &ap[kc], &c__1);
#line 257 "csptrf.f"
	    i__1 = kc + imax - 1;
#line 257 "csptrf.f"
	    colmax = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kc + 
		    imax - 1]), abs(d__2));
#line 258 "csptrf.f"
	} else {
#line 259 "csptrf.f"
	    colmax = 0.;
#line 260 "csptrf.f"
	}

#line 262 "csptrf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

#line 266 "csptrf.f"
	    if (*info == 0) {
#line 266 "csptrf.f"
		*info = k;
#line 266 "csptrf.f"
	    }
#line 268 "csptrf.f"
	    kp = k;
#line 269 "csptrf.f"
	} else {
#line 270 "csptrf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 274 "csptrf.f"
		kp = k;
#line 275 "csptrf.f"
	    } else {

#line 277 "csptrf.f"
		rowmax = 0.;
#line 278 "csptrf.f"
		jmax = imax;
#line 279 "csptrf.f"
		kx = imax * (imax + 1) / 2 + imax;
#line 280 "csptrf.f"
		i__1 = k;
#line 280 "csptrf.f"
		for (j = imax + 1; j <= i__1; ++j) {
#line 281 "csptrf.f"
		    i__2 = kx;
#line 281 "csptrf.f"
		    if ((d__1 = ap[i__2].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    kx]), abs(d__2)) > rowmax) {
#line 282 "csptrf.f"
			i__2 = kx;
#line 282 "csptrf.f"
			rowmax = (d__1 = ap[i__2].r, abs(d__1)) + (d__2 = 
				d_imag(&ap[kx]), abs(d__2));
#line 283 "csptrf.f"
			jmax = j;
#line 284 "csptrf.f"
		    }
#line 285 "csptrf.f"
		    kx += j;
#line 286 "csptrf.f"
/* L20: */
#line 286 "csptrf.f"
		}
#line 287 "csptrf.f"
		kpc = (imax - 1) * imax / 2 + 1;
#line 288 "csptrf.f"
		if (imax > 1) {
#line 289 "csptrf.f"
		    i__1 = imax - 1;
#line 289 "csptrf.f"
		    jmax = icamax_(&i__1, &ap[kpc], &c__1);
/* Computing MAX */
#line 290 "csptrf.f"
		    i__1 = kpc + jmax - 1;
#line 290 "csptrf.f"
		    d__3 = rowmax, d__4 = (d__1 = ap[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&ap[kpc + jmax - 1]), abs(d__2));
#line 290 "csptrf.f"
		    rowmax = max(d__3,d__4);
#line 291 "csptrf.f"
		}

#line 293 "csptrf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 297 "csptrf.f"
		    kp = k;
#line 298 "csptrf.f"
		} else /* if(complicated condition) */ {
#line 298 "csptrf.f"
		    i__1 = kpc + imax - 1;
#line 298 "csptrf.f"
		    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    kpc + imax - 1]), abs(d__2)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 303 "csptrf.f"
			kp = imax;
#line 304 "csptrf.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 309 "csptrf.f"
			kp = imax;
#line 310 "csptrf.f"
			kstep = 2;
#line 311 "csptrf.f"
		    }
#line 311 "csptrf.f"
		}
#line 312 "csptrf.f"
	    }

#line 314 "csptrf.f"
	    kk = k - kstep + 1;
#line 315 "csptrf.f"
	    if (kstep == 2) {
#line 315 "csptrf.f"
		knc = knc - k + 1;
#line 315 "csptrf.f"
	    }
#line 317 "csptrf.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 322 "csptrf.f"
		i__1 = kp - 1;
#line 322 "csptrf.f"
		cswap_(&i__1, &ap[knc], &c__1, &ap[kpc], &c__1);
#line 323 "csptrf.f"
		kx = kpc + kp - 1;
#line 324 "csptrf.f"
		i__1 = kk - 1;
#line 324 "csptrf.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 325 "csptrf.f"
		    kx = kx + j - 1;
#line 326 "csptrf.f"
		    i__2 = knc + j - 1;
#line 326 "csptrf.f"
		    t.r = ap[i__2].r, t.i = ap[i__2].i;
#line 327 "csptrf.f"
		    i__2 = knc + j - 1;
#line 327 "csptrf.f"
		    i__3 = kx;
#line 327 "csptrf.f"
		    ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
#line 328 "csptrf.f"
		    i__2 = kx;
#line 328 "csptrf.f"
		    ap[i__2].r = t.r, ap[i__2].i = t.i;
#line 329 "csptrf.f"
/* L30: */
#line 329 "csptrf.f"
		}
#line 330 "csptrf.f"
		i__1 = knc + kk - 1;
#line 330 "csptrf.f"
		t.r = ap[i__1].r, t.i = ap[i__1].i;
#line 331 "csptrf.f"
		i__1 = knc + kk - 1;
#line 331 "csptrf.f"
		i__2 = kpc + kp - 1;
#line 331 "csptrf.f"
		ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 332 "csptrf.f"
		i__1 = kpc + kp - 1;
#line 332 "csptrf.f"
		ap[i__1].r = t.r, ap[i__1].i = t.i;
#line 333 "csptrf.f"
		if (kstep == 2) {
#line 334 "csptrf.f"
		    i__1 = kc + k - 2;
#line 334 "csptrf.f"
		    t.r = ap[i__1].r, t.i = ap[i__1].i;
#line 335 "csptrf.f"
		    i__1 = kc + k - 2;
#line 335 "csptrf.f"
		    i__2 = kc + kp - 1;
#line 335 "csptrf.f"
		    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 336 "csptrf.f"
		    i__1 = kc + kp - 1;
#line 336 "csptrf.f"
		    ap[i__1].r = t.r, ap[i__1].i = t.i;
#line 337 "csptrf.f"
		}
#line 338 "csptrf.f"
	    }

/*           Update the leading submatrix */

#line 342 "csptrf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */

#line 354 "csptrf.f"
		z_div(&z__1, &c_b1, &ap[kc + k - 1]);
#line 354 "csptrf.f"
		r1.r = z__1.r, r1.i = z__1.i;
#line 355 "csptrf.f"
		i__1 = k - 1;
#line 355 "csptrf.f"
		z__1.r = -r1.r, z__1.i = -r1.i;
#line 355 "csptrf.f"
		cspr_(uplo, &i__1, &z__1, &ap[kc], &c__1, &ap[1], (ftnlen)1);

/*              Store U(k) in column k */

#line 359 "csptrf.f"
		i__1 = k - 1;
#line 359 "csptrf.f"
		cscal_(&i__1, &r1, &ap[kc], &c__1);
#line 360 "csptrf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */

#line 374 "csptrf.f"
		if (k > 2) {

#line 376 "csptrf.f"
		    i__1 = k - 1 + (k - 1) * k / 2;
#line 376 "csptrf.f"
		    d12.r = ap[i__1].r, d12.i = ap[i__1].i;
#line 377 "csptrf.f"
		    z_div(&z__1, &ap[k - 1 + (k - 2) * (k - 1) / 2], &d12);
#line 377 "csptrf.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 378 "csptrf.f"
		    z_div(&z__1, &ap[k + (k - 1) * k / 2], &d12);
#line 378 "csptrf.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 379 "csptrf.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 379 "csptrf.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 379 "csptrf.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 379 "csptrf.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 380 "csptrf.f"
		    z_div(&z__1, &t, &d12);
#line 380 "csptrf.f"
		    d12.r = z__1.r, d12.i = z__1.i;

#line 382 "csptrf.f"
		    for (j = k - 2; j >= 1; --j) {
#line 383 "csptrf.f"
			i__1 = j + (k - 2) * (k - 1) / 2;
#line 383 "csptrf.f"
			z__3.r = d11.r * ap[i__1].r - d11.i * ap[i__1].i, 
				z__3.i = d11.r * ap[i__1].i + d11.i * ap[i__1]
				.r;
#line 383 "csptrf.f"
			i__2 = j + (k - 1) * k / 2;
#line 383 "csptrf.f"
			z__2.r = z__3.r - ap[i__2].r, z__2.i = z__3.i - ap[
				i__2].i;
#line 383 "csptrf.f"
			z__1.r = d12.r * z__2.r - d12.i * z__2.i, z__1.i = 
				d12.r * z__2.i + d12.i * z__2.r;
#line 383 "csptrf.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 385 "csptrf.f"
			i__1 = j + (k - 1) * k / 2;
#line 385 "csptrf.f"
			z__3.r = d22.r * ap[i__1].r - d22.i * ap[i__1].i, 
				z__3.i = d22.r * ap[i__1].i + d22.i * ap[i__1]
				.r;
#line 385 "csptrf.f"
			i__2 = j + (k - 2) * (k - 1) / 2;
#line 385 "csptrf.f"
			z__2.r = z__3.r - ap[i__2].r, z__2.i = z__3.i - ap[
				i__2].i;
#line 385 "csptrf.f"
			z__1.r = d12.r * z__2.r - d12.i * z__2.i, z__1.i = 
				d12.r * z__2.i + d12.i * z__2.r;
#line 385 "csptrf.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 387 "csptrf.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 388 "csptrf.f"
			    i__1 = i__ + (j - 1) * j / 2;
#line 388 "csptrf.f"
			    i__2 = i__ + (j - 1) * j / 2;
#line 388 "csptrf.f"
			    i__3 = i__ + (k - 1) * k / 2;
#line 388 "csptrf.f"
			    z__3.r = ap[i__3].r * wk.r - ap[i__3].i * wk.i, 
				    z__3.i = ap[i__3].r * wk.i + ap[i__3].i * 
				    wk.r;
#line 388 "csptrf.f"
			    z__2.r = ap[i__2].r - z__3.r, z__2.i = ap[i__2].i 
				    - z__3.i;
#line 388 "csptrf.f"
			    i__4 = i__ + (k - 2) * (k - 1) / 2;
#line 388 "csptrf.f"
			    z__4.r = ap[i__4].r * wkm1.r - ap[i__4].i * 
				    wkm1.i, z__4.i = ap[i__4].r * wkm1.i + ap[
				    i__4].i * wkm1.r;
#line 388 "csptrf.f"
			    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - 
				    z__4.i;
#line 388 "csptrf.f"
			    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 391 "csptrf.f"
/* L40: */
#line 391 "csptrf.f"
			}
#line 392 "csptrf.f"
			i__1 = j + (k - 1) * k / 2;
#line 392 "csptrf.f"
			ap[i__1].r = wk.r, ap[i__1].i = wk.i;
#line 393 "csptrf.f"
			i__1 = j + (k - 2) * (k - 1) / 2;
#line 393 "csptrf.f"
			ap[i__1].r = wkm1.r, ap[i__1].i = wkm1.i;
#line 394 "csptrf.f"
/* L50: */
#line 394 "csptrf.f"
		    }

#line 396 "csptrf.f"
		}
#line 397 "csptrf.f"
	    }
#line 398 "csptrf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 402 "csptrf.f"
	if (kstep == 1) {
#line 403 "csptrf.f"
	    ipiv[k] = kp;
#line 404 "csptrf.f"
	} else {
#line 405 "csptrf.f"
	    ipiv[k] = -kp;
#line 406 "csptrf.f"
	    ipiv[k - 1] = -kp;
#line 407 "csptrf.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 411 "csptrf.f"
	k -= kstep;
#line 412 "csptrf.f"
	kc = knc - k;
#line 413 "csptrf.f"
	goto L10;

#line 415 "csptrf.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 422 "csptrf.f"
	k = 1;
#line 423 "csptrf.f"
	kc = 1;
#line 424 "csptrf.f"
	npp = *n * (*n + 1) / 2;
#line 425 "csptrf.f"
L60:
#line 426 "csptrf.f"
	knc = kc;

/*        If K > N, exit from loop */

#line 430 "csptrf.f"
	if (k > *n) {
#line 430 "csptrf.f"
	    goto L110;
#line 430 "csptrf.f"
	}
#line 432 "csptrf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 437 "csptrf.f"
	i__1 = kc;
#line 437 "csptrf.f"
	absakk = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kc]), 
		abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value */

#line 442 "csptrf.f"
	if (k < *n) {
#line 443 "csptrf.f"
	    i__1 = *n - k;
#line 443 "csptrf.f"
	    imax = k + icamax_(&i__1, &ap[kc + 1], &c__1);
#line 444 "csptrf.f"
	    i__1 = kc + imax - k;
#line 444 "csptrf.f"
	    colmax = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kc + 
		    imax - k]), abs(d__2));
#line 445 "csptrf.f"
	} else {
#line 446 "csptrf.f"
	    colmax = 0.;
#line 447 "csptrf.f"
	}

#line 449 "csptrf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero: set INFO and continue */

#line 453 "csptrf.f"
	    if (*info == 0) {
#line 453 "csptrf.f"
		*info = k;
#line 453 "csptrf.f"
	    }
#line 455 "csptrf.f"
	    kp = k;
#line 456 "csptrf.f"
	} else {
#line 457 "csptrf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 461 "csptrf.f"
		kp = k;
#line 462 "csptrf.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 467 "csptrf.f"
		rowmax = 0.;
#line 468 "csptrf.f"
		kx = kc + imax - k;
#line 469 "csptrf.f"
		i__1 = imax - 1;
#line 469 "csptrf.f"
		for (j = k; j <= i__1; ++j) {
#line 470 "csptrf.f"
		    i__2 = kx;
#line 470 "csptrf.f"
		    if ((d__1 = ap[i__2].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    kx]), abs(d__2)) > rowmax) {
#line 471 "csptrf.f"
			i__2 = kx;
#line 471 "csptrf.f"
			rowmax = (d__1 = ap[i__2].r, abs(d__1)) + (d__2 = 
				d_imag(&ap[kx]), abs(d__2));
#line 472 "csptrf.f"
			jmax = j;
#line 473 "csptrf.f"
		    }
#line 474 "csptrf.f"
		    kx = kx + *n - j;
#line 475 "csptrf.f"
/* L70: */
#line 475 "csptrf.f"
		}
#line 476 "csptrf.f"
		kpc = npp - (*n - imax + 1) * (*n - imax + 2) / 2 + 1;
#line 477 "csptrf.f"
		if (imax < *n) {
#line 478 "csptrf.f"
		    i__1 = *n - imax;
#line 478 "csptrf.f"
		    jmax = imax + icamax_(&i__1, &ap[kpc + 1], &c__1);
/* Computing MAX */
#line 479 "csptrf.f"
		    i__1 = kpc + jmax - imax;
#line 479 "csptrf.f"
		    d__3 = rowmax, d__4 = (d__1 = ap[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&ap[kpc + jmax - imax]), abs(d__2));
#line 479 "csptrf.f"
		    rowmax = max(d__3,d__4);
#line 480 "csptrf.f"
		}

#line 482 "csptrf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 486 "csptrf.f"
		    kp = k;
#line 487 "csptrf.f"
		} else /* if(complicated condition) */ {
#line 487 "csptrf.f"
		    i__1 = kpc;
#line 487 "csptrf.f"
		    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    kpc]), abs(d__2)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 492 "csptrf.f"
			kp = imax;
#line 493 "csptrf.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 498 "csptrf.f"
			kp = imax;
#line 499 "csptrf.f"
			kstep = 2;
#line 500 "csptrf.f"
		    }
#line 500 "csptrf.f"
		}
#line 501 "csptrf.f"
	    }

#line 503 "csptrf.f"
	    kk = k + kstep - 1;
#line 504 "csptrf.f"
	    if (kstep == 2) {
#line 504 "csptrf.f"
		knc = knc + *n - k + 1;
#line 504 "csptrf.f"
	    }
#line 506 "csptrf.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 511 "csptrf.f"
		if (kp < *n) {
#line 511 "csptrf.f"
		    i__1 = *n - kp;
#line 511 "csptrf.f"
		    cswap_(&i__1, &ap[knc + kp - kk + 1], &c__1, &ap[kpc + 1],
			     &c__1);
#line 511 "csptrf.f"
		}
#line 514 "csptrf.f"
		kx = knc + kp - kk;
#line 515 "csptrf.f"
		i__1 = kp - 1;
#line 515 "csptrf.f"
		for (j = kk + 1; j <= i__1; ++j) {
#line 516 "csptrf.f"
		    kx = kx + *n - j + 1;
#line 517 "csptrf.f"
		    i__2 = knc + j - kk;
#line 517 "csptrf.f"
		    t.r = ap[i__2].r, t.i = ap[i__2].i;
#line 518 "csptrf.f"
		    i__2 = knc + j - kk;
#line 518 "csptrf.f"
		    i__3 = kx;
#line 518 "csptrf.f"
		    ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
#line 519 "csptrf.f"
		    i__2 = kx;
#line 519 "csptrf.f"
		    ap[i__2].r = t.r, ap[i__2].i = t.i;
#line 520 "csptrf.f"
/* L80: */
#line 520 "csptrf.f"
		}
#line 521 "csptrf.f"
		i__1 = knc;
#line 521 "csptrf.f"
		t.r = ap[i__1].r, t.i = ap[i__1].i;
#line 522 "csptrf.f"
		i__1 = knc;
#line 522 "csptrf.f"
		i__2 = kpc;
#line 522 "csptrf.f"
		ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 523 "csptrf.f"
		i__1 = kpc;
#line 523 "csptrf.f"
		ap[i__1].r = t.r, ap[i__1].i = t.i;
#line 524 "csptrf.f"
		if (kstep == 2) {
#line 525 "csptrf.f"
		    i__1 = kc + 1;
#line 525 "csptrf.f"
		    t.r = ap[i__1].r, t.i = ap[i__1].i;
#line 526 "csptrf.f"
		    i__1 = kc + 1;
#line 526 "csptrf.f"
		    i__2 = kc + kp - k;
#line 526 "csptrf.f"
		    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 527 "csptrf.f"
		    i__1 = kc + kp - k;
#line 527 "csptrf.f"
		    ap[i__1].r = t.r, ap[i__1].i = t.i;
#line 528 "csptrf.f"
		}
#line 529 "csptrf.f"
	    }

/*           Update the trailing submatrix */

#line 533 "csptrf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 541 "csptrf.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */

#line 547 "csptrf.f"
		    z_div(&z__1, &c_b1, &ap[kc]);
#line 547 "csptrf.f"
		    r1.r = z__1.r, r1.i = z__1.i;
#line 548 "csptrf.f"
		    i__1 = *n - k;
#line 548 "csptrf.f"
		    z__1.r = -r1.r, z__1.i = -r1.i;
#line 548 "csptrf.f"
		    cspr_(uplo, &i__1, &z__1, &ap[kc + 1], &c__1, &ap[kc + *n 
			    - k + 1], (ftnlen)1);

/*                 Store L(k) in column K */

#line 553 "csptrf.f"
		    i__1 = *n - k;
#line 553 "csptrf.f"
		    cscal_(&i__1, &r1, &ap[kc + 1], &c__1);
#line 554 "csptrf.f"
		}
#line 555 "csptrf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns K and K+1 now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

#line 564 "csptrf.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T */
/*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 574 "csptrf.f"
		    i__1 = k + 1 + (k - 1) * ((*n << 1) - k) / 2;
#line 574 "csptrf.f"
		    d21.r = ap[i__1].r, d21.i = ap[i__1].i;
#line 575 "csptrf.f"
		    z_div(&z__1, &ap[k + 1 + k * ((*n << 1) - k - 1) / 2], &
			    d21);
#line 575 "csptrf.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 576 "csptrf.f"
		    z_div(&z__1, &ap[k + (k - 1) * ((*n << 1) - k) / 2], &d21)
			    ;
#line 576 "csptrf.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 577 "csptrf.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 577 "csptrf.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 577 "csptrf.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 577 "csptrf.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 578 "csptrf.f"
		    z_div(&z__1, &t, &d21);
#line 578 "csptrf.f"
		    d21.r = z__1.r, d21.i = z__1.i;

#line 580 "csptrf.f"
		    i__1 = *n;
#line 580 "csptrf.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 581 "csptrf.f"
			i__2 = j + (k - 1) * ((*n << 1) - k) / 2;
#line 581 "csptrf.f"
			z__3.r = d11.r * ap[i__2].r - d11.i * ap[i__2].i, 
				z__3.i = d11.r * ap[i__2].i + d11.i * ap[i__2]
				.r;
#line 581 "csptrf.f"
			i__3 = j + k * ((*n << 1) - k - 1) / 2;
#line 581 "csptrf.f"
			z__2.r = z__3.r - ap[i__3].r, z__2.i = z__3.i - ap[
				i__3].i;
#line 581 "csptrf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 581 "csptrf.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 583 "csptrf.f"
			i__2 = j + k * ((*n << 1) - k - 1) / 2;
#line 583 "csptrf.f"
			z__3.r = d22.r * ap[i__2].r - d22.i * ap[i__2].i, 
				z__3.i = d22.r * ap[i__2].i + d22.i * ap[i__2]
				.r;
#line 583 "csptrf.f"
			i__3 = j + (k - 1) * ((*n << 1) - k) / 2;
#line 583 "csptrf.f"
			z__2.r = z__3.r - ap[i__3].r, z__2.i = z__3.i - ap[
				i__3].i;
#line 583 "csptrf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 583 "csptrf.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;
#line 585 "csptrf.f"
			i__2 = *n;
#line 585 "csptrf.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 586 "csptrf.f"
			    i__3 = i__ + (j - 1) * ((*n << 1) - j) / 2;
#line 586 "csptrf.f"
			    i__4 = i__ + (j - 1) * ((*n << 1) - j) / 2;
#line 586 "csptrf.f"
			    i__5 = i__ + (k - 1) * ((*n << 1) - k) / 2;
#line 586 "csptrf.f"
			    z__3.r = ap[i__5].r * wk.r - ap[i__5].i * wk.i, 
				    z__3.i = ap[i__5].r * wk.i + ap[i__5].i * 
				    wk.r;
#line 586 "csptrf.f"
			    z__2.r = ap[i__4].r - z__3.r, z__2.i = ap[i__4].i 
				    - z__3.i;
#line 586 "csptrf.f"
			    i__6 = i__ + k * ((*n << 1) - k - 1) / 2;
#line 586 "csptrf.f"
			    z__4.r = ap[i__6].r * wkp1.r - ap[i__6].i * 
				    wkp1.i, z__4.i = ap[i__6].r * wkp1.i + ap[
				    i__6].i * wkp1.r;
#line 586 "csptrf.f"
			    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - 
				    z__4.i;
#line 586 "csptrf.f"
			    ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 589 "csptrf.f"
/* L90: */
#line 589 "csptrf.f"
			}
#line 590 "csptrf.f"
			i__2 = j + (k - 1) * ((*n << 1) - k) / 2;
#line 590 "csptrf.f"
			ap[i__2].r = wk.r, ap[i__2].i = wk.i;
#line 591 "csptrf.f"
			i__2 = j + k * ((*n << 1) - k - 1) / 2;
#line 591 "csptrf.f"
			ap[i__2].r = wkp1.r, ap[i__2].i = wkp1.i;
#line 592 "csptrf.f"
/* L100: */
#line 592 "csptrf.f"
		    }
#line 593 "csptrf.f"
		}
#line 594 "csptrf.f"
	    }
#line 595 "csptrf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 599 "csptrf.f"
	if (kstep == 1) {
#line 600 "csptrf.f"
	    ipiv[k] = kp;
#line 601 "csptrf.f"
	} else {
#line 602 "csptrf.f"
	    ipiv[k] = -kp;
#line 603 "csptrf.f"
	    ipiv[k + 1] = -kp;
#line 604 "csptrf.f"
	}

/*        Increase K and return to the start of the main loop */

#line 608 "csptrf.f"
	k += kstep;
#line 609 "csptrf.f"
	kc = knc + *n - k + 2;
#line 610 "csptrf.f"
	goto L60;

#line 612 "csptrf.f"
    }

#line 614 "csptrf.f"
L110:
#line 615 "csptrf.f"
    return 0;

/*     End of CSPTRF */

} /* csptrf_ */

