#line 1 "zsytf2.f"
/* zsytf2.f -- translated by f2c (version 20100827).
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

#line 1 "zsytf2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal piv
oting method (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTF2( UPLO, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYTF2 computes the factorization of a complex symmetric matrix A */
/* > using the Bunch-Kaufman diagonal pivoting method: */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \ingroup complex16SYcomputational */

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
/* >    Replace l.209 and l.377 */
/* >         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN */
/* >    by */
/* >         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN */
/* > */
/* >  1-96 - Based on modifications by J. Lewis, Boeing Computer Services */
/* >         Company */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zsytf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex t, r1, d11, d12, d21, d22;
    static integer kk, kp;
    static doublecomplex wk, wkm1, wkp1;
    static integer imax, jmax;
    extern /* Subroutine */ int zsyr_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal absakk;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal colmax;
    extern integer izamax_(integer *, doublecomplex *, integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 245 "zsytf2.f"
    /* Parameter adjustments */
#line 245 "zsytf2.f"
    a_dim1 = *lda;
#line 245 "zsytf2.f"
    a_offset = 1 + a_dim1;
#line 245 "zsytf2.f"
    a -= a_offset;
#line 245 "zsytf2.f"
    --ipiv;
#line 245 "zsytf2.f"

#line 245 "zsytf2.f"
    /* Function Body */
#line 245 "zsytf2.f"
    *info = 0;
#line 246 "zsytf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 247 "zsytf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 248 "zsytf2.f"
	*info = -1;
#line 249 "zsytf2.f"
    } else if (*n < 0) {
#line 250 "zsytf2.f"
	*info = -2;
#line 251 "zsytf2.f"
    } else if (*lda < max(1,*n)) {
#line 252 "zsytf2.f"
	*info = -4;
#line 253 "zsytf2.f"
    }
#line 254 "zsytf2.f"
    if (*info != 0) {
#line 255 "zsytf2.f"
	i__1 = -(*info);
#line 255 "zsytf2.f"
	xerbla_("ZSYTF2", &i__1, (ftnlen)6);
#line 256 "zsytf2.f"
	return 0;
#line 257 "zsytf2.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 261 "zsytf2.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 263 "zsytf2.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 270 "zsytf2.f"
	k = *n;
#line 271 "zsytf2.f"
L10:

/*        If K < 1, exit from loop */

#line 275 "zsytf2.f"
	if (k < 1) {
#line 275 "zsytf2.f"
	    goto L70;
#line 275 "zsytf2.f"
	}
#line 277 "zsytf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 282 "zsytf2.f"
	i__1 = k + k * a_dim1;
#line 282 "zsytf2.f"
	absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * 
		a_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 288 "zsytf2.f"
	if (k > 1) {
#line 289 "zsytf2.f"
	    i__1 = k - 1;
#line 289 "zsytf2.f"
	    imax = izamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 290 "zsytf2.f"
	    i__1 = imax + k * a_dim1;
#line 290 "zsytf2.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 291 "zsytf2.f"
	} else {
#line 292 "zsytf2.f"
	    colmax = 0.;
#line 293 "zsytf2.f"
	}

#line 295 "zsytf2.f"
	if (max(absakk,colmax) == 0. || disnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 300 "zsytf2.f"
	    if (*info == 0) {
#line 300 "zsytf2.f"
		*info = k;
#line 300 "zsytf2.f"
	    }
#line 302 "zsytf2.f"
	    kp = k;
#line 303 "zsytf2.f"
	} else {
#line 304 "zsytf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 308 "zsytf2.f"
		kp = k;
#line 309 "zsytf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 314 "zsytf2.f"
		i__1 = k - imax;
#line 314 "zsytf2.f"
		jmax = imax + izamax_(&i__1, &a[imax + (imax + 1) * a_dim1], 
			lda);
#line 315 "zsytf2.f"
		i__1 = imax + jmax * a_dim1;
#line 315 "zsytf2.f"
		rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			imax + jmax * a_dim1]), abs(d__2));
#line 316 "zsytf2.f"
		if (imax > 1) {
#line 317 "zsytf2.f"
		    i__1 = imax - 1;
#line 317 "zsytf2.f"
		    jmax = izamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
#line 318 "zsytf2.f"
		    i__1 = jmax + imax * a_dim1;
#line 318 "zsytf2.f"
		    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&a[jmax + imax * a_dim1]), abs(d__2)
			    );
#line 318 "zsytf2.f"
		    rowmax = max(d__3,d__4);
#line 319 "zsytf2.f"
		}

#line 321 "zsytf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 325 "zsytf2.f"
		    kp = k;
#line 326 "zsytf2.f"
		} else /* if(complicated condition) */ {
#line 326 "zsytf2.f"
		    i__1 = imax + imax * a_dim1;
#line 326 "zsytf2.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    imax + imax * a_dim1]), abs(d__2)) >= alpha * 
			    rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 331 "zsytf2.f"
			kp = imax;
#line 332 "zsytf2.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 337 "zsytf2.f"
			kp = imax;
#line 338 "zsytf2.f"
			kstep = 2;
#line 339 "zsytf2.f"
		    }
#line 339 "zsytf2.f"
		}
#line 340 "zsytf2.f"
	    }

#line 342 "zsytf2.f"
	    kk = k - kstep + 1;
#line 343 "zsytf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 348 "zsytf2.f"
		i__1 = kp - 1;
#line 348 "zsytf2.f"
		zswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1],
			 &c__1);
#line 349 "zsytf2.f"
		i__1 = kk - kp - 1;
#line 349 "zsytf2.f"
		zswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 351 "zsytf2.f"
		i__1 = kk + kk * a_dim1;
#line 351 "zsytf2.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 352 "zsytf2.f"
		i__1 = kk + kk * a_dim1;
#line 352 "zsytf2.f"
		i__2 = kp + kp * a_dim1;
#line 352 "zsytf2.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 353 "zsytf2.f"
		i__1 = kp + kp * a_dim1;
#line 353 "zsytf2.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 354 "zsytf2.f"
		if (kstep == 2) {
#line 355 "zsytf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 355 "zsytf2.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 356 "zsytf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 356 "zsytf2.f"
		    i__2 = kp + k * a_dim1;
#line 356 "zsytf2.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 357 "zsytf2.f"
		    i__1 = kp + k * a_dim1;
#line 357 "zsytf2.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 358 "zsytf2.f"
		}
#line 359 "zsytf2.f"
	    }

/*           Update the leading submatrix */

#line 363 "zsytf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */

#line 375 "zsytf2.f"
		z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 375 "zsytf2.f"
		r1.r = z__1.r, r1.i = z__1.i;
#line 376 "zsytf2.f"
		i__1 = k - 1;
#line 376 "zsytf2.f"
		z__1.r = -r1.r, z__1.i = -r1.i;
#line 376 "zsytf2.f"
		zsyr_(uplo, &i__1, &z__1, &a[k * a_dim1 + 1], &c__1, &a[
			a_offset], lda, (ftnlen)1);

/*              Store U(k) in column k */

#line 380 "zsytf2.f"
		i__1 = k - 1;
#line 380 "zsytf2.f"
		zscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 381 "zsytf2.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */

#line 395 "zsytf2.f"
		if (k > 2) {

#line 397 "zsytf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 397 "zsytf2.f"
		    d12.r = a[i__1].r, d12.i = a[i__1].i;
#line 398 "zsytf2.f"
		    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &d12);
#line 398 "zsytf2.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 399 "zsytf2.f"
		    z_div(&z__1, &a[k + k * a_dim1], &d12);
#line 399 "zsytf2.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 400 "zsytf2.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 400 "zsytf2.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 400 "zsytf2.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 400 "zsytf2.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 401 "zsytf2.f"
		    z_div(&z__1, &t, &d12);
#line 401 "zsytf2.f"
		    d12.r = z__1.r, d12.i = z__1.i;

#line 403 "zsytf2.f"
		    for (j = k - 2; j >= 1; --j) {
#line 404 "zsytf2.f"
			i__1 = j + (k - 1) * a_dim1;
#line 404 "zsytf2.f"
			z__3.r = d11.r * a[i__1].r - d11.i * a[i__1].i, 
				z__3.i = d11.r * a[i__1].i + d11.i * a[i__1]
				.r;
#line 404 "zsytf2.f"
			i__2 = j + k * a_dim1;
#line 404 "zsytf2.f"
			z__2.r = z__3.r - a[i__2].r, z__2.i = z__3.i - a[i__2]
				.i;
#line 404 "zsytf2.f"
			z__1.r = d12.r * z__2.r - d12.i * z__2.i, z__1.i = 
				d12.r * z__2.i + d12.i * z__2.r;
#line 404 "zsytf2.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 405 "zsytf2.f"
			i__1 = j + k * a_dim1;
#line 405 "zsytf2.f"
			z__3.r = d22.r * a[i__1].r - d22.i * a[i__1].i, 
				z__3.i = d22.r * a[i__1].i + d22.i * a[i__1]
				.r;
#line 405 "zsytf2.f"
			i__2 = j + (k - 1) * a_dim1;
#line 405 "zsytf2.f"
			z__2.r = z__3.r - a[i__2].r, z__2.i = z__3.i - a[i__2]
				.i;
#line 405 "zsytf2.f"
			z__1.r = d12.r * z__2.r - d12.i * z__2.i, z__1.i = 
				d12.r * z__2.i + d12.i * z__2.r;
#line 405 "zsytf2.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 406 "zsytf2.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 407 "zsytf2.f"
			    i__1 = i__ + j * a_dim1;
#line 407 "zsytf2.f"
			    i__2 = i__ + j * a_dim1;
#line 407 "zsytf2.f"
			    i__3 = i__ + k * a_dim1;
#line 407 "zsytf2.f"
			    z__3.r = a[i__3].r * wk.r - a[i__3].i * wk.i, 
				    z__3.i = a[i__3].r * wk.i + a[i__3].i * 
				    wk.r;
#line 407 "zsytf2.f"
			    z__2.r = a[i__2].r - z__3.r, z__2.i = a[i__2].i - 
				    z__3.i;
#line 407 "zsytf2.f"
			    i__4 = i__ + (k - 1) * a_dim1;
#line 407 "zsytf2.f"
			    z__4.r = a[i__4].r * wkm1.r - a[i__4].i * wkm1.i, 
				    z__4.i = a[i__4].r * wkm1.i + a[i__4].i * 
				    wkm1.r;
#line 407 "zsytf2.f"
			    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - 
				    z__4.i;
#line 407 "zsytf2.f"
			    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 409 "zsytf2.f"
/* L20: */
#line 409 "zsytf2.f"
			}
#line 410 "zsytf2.f"
			i__1 = j + k * a_dim1;
#line 410 "zsytf2.f"
			a[i__1].r = wk.r, a[i__1].i = wk.i;
#line 411 "zsytf2.f"
			i__1 = j + (k - 1) * a_dim1;
#line 411 "zsytf2.f"
			a[i__1].r = wkm1.r, a[i__1].i = wkm1.i;
#line 412 "zsytf2.f"
/* L30: */
#line 412 "zsytf2.f"
		    }

#line 414 "zsytf2.f"
		}

#line 416 "zsytf2.f"
	    }
#line 417 "zsytf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 421 "zsytf2.f"
	if (kstep == 1) {
#line 422 "zsytf2.f"
	    ipiv[k] = kp;
#line 423 "zsytf2.f"
	} else {
#line 424 "zsytf2.f"
	    ipiv[k] = -kp;
#line 425 "zsytf2.f"
	    ipiv[k - 1] = -kp;
#line 426 "zsytf2.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 430 "zsytf2.f"
	k -= kstep;
#line 431 "zsytf2.f"
	goto L10;

#line 433 "zsytf2.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 440 "zsytf2.f"
	k = 1;
#line 441 "zsytf2.f"
L40:

/*        If K > N, exit from loop */

#line 445 "zsytf2.f"
	if (k > *n) {
#line 445 "zsytf2.f"
	    goto L70;
#line 445 "zsytf2.f"
	}
#line 447 "zsytf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 452 "zsytf2.f"
	i__1 = k + k * a_dim1;
#line 452 "zsytf2.f"
	absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * 
		a_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 458 "zsytf2.f"
	if (k < *n) {
#line 459 "zsytf2.f"
	    i__1 = *n - k;
#line 459 "zsytf2.f"
	    imax = k + izamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 460 "zsytf2.f"
	    i__1 = imax + k * a_dim1;
#line 460 "zsytf2.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 461 "zsytf2.f"
	} else {
#line 462 "zsytf2.f"
	    colmax = 0.;
#line 463 "zsytf2.f"
	}

#line 465 "zsytf2.f"
	if (max(absakk,colmax) == 0. || disnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 470 "zsytf2.f"
	    if (*info == 0) {
#line 470 "zsytf2.f"
		*info = k;
#line 470 "zsytf2.f"
	    }
#line 472 "zsytf2.f"
	    kp = k;
#line 473 "zsytf2.f"
	} else {
#line 474 "zsytf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 478 "zsytf2.f"
		kp = k;
#line 479 "zsytf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 484 "zsytf2.f"
		i__1 = imax - k;
#line 484 "zsytf2.f"
		jmax = k - 1 + izamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 485 "zsytf2.f"
		i__1 = imax + jmax * a_dim1;
#line 485 "zsytf2.f"
		rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			imax + jmax * a_dim1]), abs(d__2));
#line 486 "zsytf2.f"
		if (imax < *n) {
#line 487 "zsytf2.f"
		    i__1 = *n - imax;
#line 487 "zsytf2.f"
		    jmax = imax + izamax_(&i__1, &a[imax + 1 + imax * a_dim1],
			     &c__1);
/* Computing MAX */
#line 488 "zsytf2.f"
		    i__1 = jmax + imax * a_dim1;
#line 488 "zsytf2.f"
		    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&a[jmax + imax * a_dim1]), abs(d__2)
			    );
#line 488 "zsytf2.f"
		    rowmax = max(d__3,d__4);
#line 489 "zsytf2.f"
		}

#line 491 "zsytf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 495 "zsytf2.f"
		    kp = k;
#line 496 "zsytf2.f"
		} else /* if(complicated condition) */ {
#line 496 "zsytf2.f"
		    i__1 = imax + imax * a_dim1;
#line 496 "zsytf2.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    imax + imax * a_dim1]), abs(d__2)) >= alpha * 
			    rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 501 "zsytf2.f"
			kp = imax;
#line 502 "zsytf2.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 507 "zsytf2.f"
			kp = imax;
#line 508 "zsytf2.f"
			kstep = 2;
#line 509 "zsytf2.f"
		    }
#line 509 "zsytf2.f"
		}
#line 510 "zsytf2.f"
	    }

#line 512 "zsytf2.f"
	    kk = k + kstep - 1;
#line 513 "zsytf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 518 "zsytf2.f"
		if (kp < *n) {
#line 518 "zsytf2.f"
		    i__1 = *n - kp;
#line 518 "zsytf2.f"
		    zswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 518 "zsytf2.f"
		}
#line 520 "zsytf2.f"
		i__1 = kp - kk - 1;
#line 520 "zsytf2.f"
		zswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 522 "zsytf2.f"
		i__1 = kk + kk * a_dim1;
#line 522 "zsytf2.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 523 "zsytf2.f"
		i__1 = kk + kk * a_dim1;
#line 523 "zsytf2.f"
		i__2 = kp + kp * a_dim1;
#line 523 "zsytf2.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 524 "zsytf2.f"
		i__1 = kp + kp * a_dim1;
#line 524 "zsytf2.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 525 "zsytf2.f"
		if (kstep == 2) {
#line 526 "zsytf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 526 "zsytf2.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 527 "zsytf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 527 "zsytf2.f"
		    i__2 = kp + k * a_dim1;
#line 527 "zsytf2.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 528 "zsytf2.f"
		    i__1 = kp + k * a_dim1;
#line 528 "zsytf2.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 529 "zsytf2.f"
		}
#line 530 "zsytf2.f"
	    }

/*           Update the trailing submatrix */

#line 534 "zsytf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 542 "zsytf2.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */

#line 548 "zsytf2.f"
		    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 548 "zsytf2.f"
		    r1.r = z__1.r, r1.i = z__1.i;
#line 549 "zsytf2.f"
		    i__1 = *n - k;
#line 549 "zsytf2.f"
		    z__1.r = -r1.r, z__1.i = -r1.i;
#line 549 "zsytf2.f"
		    zsyr_(uplo, &i__1, &z__1, &a[k + 1 + k * a_dim1], &c__1, &
			    a[k + 1 + (k + 1) * a_dim1], lda, (ftnlen)1);

/*                 Store L(k) in column K */

#line 554 "zsytf2.f"
		    i__1 = *n - k;
#line 554 "zsytf2.f"
		    zscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 555 "zsytf2.f"
		}
#line 556 "zsytf2.f"
	    } else {

/*              2-by-2 pivot block D(k) */

#line 560 "zsytf2.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T */
/*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 570 "zsytf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 570 "zsytf2.f"
		    d21.r = a[i__1].r, d21.i = a[i__1].i;
#line 571 "zsytf2.f"
		    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &d21);
#line 571 "zsytf2.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 572 "zsytf2.f"
		    z_div(&z__1, &a[k + k * a_dim1], &d21);
#line 572 "zsytf2.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 573 "zsytf2.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 573 "zsytf2.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 573 "zsytf2.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 573 "zsytf2.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 574 "zsytf2.f"
		    z_div(&z__1, &t, &d21);
#line 574 "zsytf2.f"
		    d21.r = z__1.r, d21.i = z__1.i;

#line 576 "zsytf2.f"
		    i__1 = *n;
#line 576 "zsytf2.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 577 "zsytf2.f"
			i__2 = j + k * a_dim1;
#line 577 "zsytf2.f"
			z__3.r = d11.r * a[i__2].r - d11.i * a[i__2].i, 
				z__3.i = d11.r * a[i__2].i + d11.i * a[i__2]
				.r;
#line 577 "zsytf2.f"
			i__3 = j + (k + 1) * a_dim1;
#line 577 "zsytf2.f"
			z__2.r = z__3.r - a[i__3].r, z__2.i = z__3.i - a[i__3]
				.i;
#line 577 "zsytf2.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 577 "zsytf2.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 578 "zsytf2.f"
			i__2 = j + (k + 1) * a_dim1;
#line 578 "zsytf2.f"
			z__3.r = d22.r * a[i__2].r - d22.i * a[i__2].i, 
				z__3.i = d22.r * a[i__2].i + d22.i * a[i__2]
				.r;
#line 578 "zsytf2.f"
			i__3 = j + k * a_dim1;
#line 578 "zsytf2.f"
			z__2.r = z__3.r - a[i__3].r, z__2.i = z__3.i - a[i__3]
				.i;
#line 578 "zsytf2.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 578 "zsytf2.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;
#line 579 "zsytf2.f"
			i__2 = *n;
#line 579 "zsytf2.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 580 "zsytf2.f"
			    i__3 = i__ + j * a_dim1;
#line 580 "zsytf2.f"
			    i__4 = i__ + j * a_dim1;
#line 580 "zsytf2.f"
			    i__5 = i__ + k * a_dim1;
#line 580 "zsytf2.f"
			    z__3.r = a[i__5].r * wk.r - a[i__5].i * wk.i, 
				    z__3.i = a[i__5].r * wk.i + a[i__5].i * 
				    wk.r;
#line 580 "zsytf2.f"
			    z__2.r = a[i__4].r - z__3.r, z__2.i = a[i__4].i - 
				    z__3.i;
#line 580 "zsytf2.f"
			    i__6 = i__ + (k + 1) * a_dim1;
#line 580 "zsytf2.f"
			    z__4.r = a[i__6].r * wkp1.r - a[i__6].i * wkp1.i, 
				    z__4.i = a[i__6].r * wkp1.i + a[i__6].i * 
				    wkp1.r;
#line 580 "zsytf2.f"
			    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - 
				    z__4.i;
#line 580 "zsytf2.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 582 "zsytf2.f"
/* L50: */
#line 582 "zsytf2.f"
			}
#line 583 "zsytf2.f"
			i__2 = j + k * a_dim1;
#line 583 "zsytf2.f"
			a[i__2].r = wk.r, a[i__2].i = wk.i;
#line 584 "zsytf2.f"
			i__2 = j + (k + 1) * a_dim1;
#line 584 "zsytf2.f"
			a[i__2].r = wkp1.r, a[i__2].i = wkp1.i;
#line 585 "zsytf2.f"
/* L60: */
#line 585 "zsytf2.f"
		    }
#line 586 "zsytf2.f"
		}
#line 587 "zsytf2.f"
	    }
#line 588 "zsytf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 592 "zsytf2.f"
	if (kstep == 1) {
#line 593 "zsytf2.f"
	    ipiv[k] = kp;
#line 594 "zsytf2.f"
	} else {
#line 595 "zsytf2.f"
	    ipiv[k] = -kp;
#line 596 "zsytf2.f"
	    ipiv[k + 1] = -kp;
#line 597 "zsytf2.f"
	}

/*        Increase K and return to the start of the main loop */

#line 601 "zsytf2.f"
	k += kstep;
#line 602 "zsytf2.f"
	goto L40;

#line 604 "zsytf2.f"
    }

#line 606 "zsytf2.f"
L70:
#line 607 "zsytf2.f"
    return 0;

/*     End of ZSYTF2 */

} /* zsytf2_ */

