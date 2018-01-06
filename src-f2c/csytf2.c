#line 1 "csytf2.f"
/* csytf2.f -- translated by f2c (version 20100827).
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

#line 1 "csytf2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal piv
oting method (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTF2( UPLO, N, A, LDA, IPIV, INFO ) */

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
/* > CSYTF2 computes the factorization of a complex symmetric matrix A */
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
/* >  09-29-06 - patch from */
/* >    Bobby Cheng, MathWorks */
/* > */
/* >    Replace l.209 and l.377 */
/* >         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN */
/* >    by */
/* >         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN */
/* > */
/* >  1-96 - Based on modifications by J. Lewis, Boeing Computer Services */
/* >         Company */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int csytf2_(char *uplo, integer *n, doublecomplex *a, 
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
    extern /* Subroutine */ int csyr_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
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
    static doublereal colmax;
    extern logical sisnan_(doublereal *);
    static doublereal rowmax;


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

#line 245 "csytf2.f"
    /* Parameter adjustments */
#line 245 "csytf2.f"
    a_dim1 = *lda;
#line 245 "csytf2.f"
    a_offset = 1 + a_dim1;
#line 245 "csytf2.f"
    a -= a_offset;
#line 245 "csytf2.f"
    --ipiv;
#line 245 "csytf2.f"

#line 245 "csytf2.f"
    /* Function Body */
#line 245 "csytf2.f"
    *info = 0;
#line 246 "csytf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 247 "csytf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 248 "csytf2.f"
	*info = -1;
#line 249 "csytf2.f"
    } else if (*n < 0) {
#line 250 "csytf2.f"
	*info = -2;
#line 251 "csytf2.f"
    } else if (*lda < max(1,*n)) {
#line 252 "csytf2.f"
	*info = -4;
#line 253 "csytf2.f"
    }
#line 254 "csytf2.f"
    if (*info != 0) {
#line 255 "csytf2.f"
	i__1 = -(*info);
#line 255 "csytf2.f"
	xerbla_("CSYTF2", &i__1, (ftnlen)6);
#line 256 "csytf2.f"
	return 0;
#line 257 "csytf2.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 261 "csytf2.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 263 "csytf2.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 270 "csytf2.f"
	k = *n;
#line 271 "csytf2.f"
L10:

/*        If K < 1, exit from loop */

#line 275 "csytf2.f"
	if (k < 1) {
#line 275 "csytf2.f"
	    goto L70;
#line 275 "csytf2.f"
	}
#line 277 "csytf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 282 "csytf2.f"
	i__1 = k + k * a_dim1;
#line 282 "csytf2.f"
	absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * 
		a_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 288 "csytf2.f"
	if (k > 1) {
#line 289 "csytf2.f"
	    i__1 = k - 1;
#line 289 "csytf2.f"
	    imax = icamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 290 "csytf2.f"
	    i__1 = imax + k * a_dim1;
#line 290 "csytf2.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 291 "csytf2.f"
	} else {
#line 292 "csytf2.f"
	    colmax = 0.;
#line 293 "csytf2.f"
	}

#line 295 "csytf2.f"
	if (max(absakk,colmax) == 0. || sisnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 300 "csytf2.f"
	    if (*info == 0) {
#line 300 "csytf2.f"
		*info = k;
#line 300 "csytf2.f"
	    }
#line 302 "csytf2.f"
	    kp = k;
#line 303 "csytf2.f"
	} else {
#line 304 "csytf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 308 "csytf2.f"
		kp = k;
#line 309 "csytf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 314 "csytf2.f"
		i__1 = k - imax;
#line 314 "csytf2.f"
		jmax = imax + icamax_(&i__1, &a[imax + (imax + 1) * a_dim1], 
			lda);
#line 315 "csytf2.f"
		i__1 = imax + jmax * a_dim1;
#line 315 "csytf2.f"
		rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			imax + jmax * a_dim1]), abs(d__2));
#line 316 "csytf2.f"
		if (imax > 1) {
#line 317 "csytf2.f"
		    i__1 = imax - 1;
#line 317 "csytf2.f"
		    jmax = icamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
#line 318 "csytf2.f"
		    i__1 = jmax + imax * a_dim1;
#line 318 "csytf2.f"
		    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&a[jmax + imax * a_dim1]), abs(d__2)
			    );
#line 318 "csytf2.f"
		    rowmax = max(d__3,d__4);
#line 319 "csytf2.f"
		}

#line 321 "csytf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 325 "csytf2.f"
		    kp = k;
#line 326 "csytf2.f"
		} else /* if(complicated condition) */ {
#line 326 "csytf2.f"
		    i__1 = imax + imax * a_dim1;
#line 326 "csytf2.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    imax + imax * a_dim1]), abs(d__2)) >= alpha * 
			    rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 331 "csytf2.f"
			kp = imax;
#line 332 "csytf2.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 337 "csytf2.f"
			kp = imax;
#line 338 "csytf2.f"
			kstep = 2;
#line 339 "csytf2.f"
		    }
#line 339 "csytf2.f"
		}
#line 340 "csytf2.f"
	    }

#line 342 "csytf2.f"
	    kk = k - kstep + 1;
#line 343 "csytf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 348 "csytf2.f"
		i__1 = kp - 1;
#line 348 "csytf2.f"
		cswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1],
			 &c__1);
#line 349 "csytf2.f"
		i__1 = kk - kp - 1;
#line 349 "csytf2.f"
		cswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 351 "csytf2.f"
		i__1 = kk + kk * a_dim1;
#line 351 "csytf2.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 352 "csytf2.f"
		i__1 = kk + kk * a_dim1;
#line 352 "csytf2.f"
		i__2 = kp + kp * a_dim1;
#line 352 "csytf2.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 353 "csytf2.f"
		i__1 = kp + kp * a_dim1;
#line 353 "csytf2.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 354 "csytf2.f"
		if (kstep == 2) {
#line 355 "csytf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 355 "csytf2.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 356 "csytf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 356 "csytf2.f"
		    i__2 = kp + k * a_dim1;
#line 356 "csytf2.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 357 "csytf2.f"
		    i__1 = kp + k * a_dim1;
#line 357 "csytf2.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 358 "csytf2.f"
		}
#line 359 "csytf2.f"
	    }

/*           Update the leading submatrix */

#line 363 "csytf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */

#line 375 "csytf2.f"
		z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 375 "csytf2.f"
		r1.r = z__1.r, r1.i = z__1.i;
#line 376 "csytf2.f"
		i__1 = k - 1;
#line 376 "csytf2.f"
		z__1.r = -r1.r, z__1.i = -r1.i;
#line 376 "csytf2.f"
		csyr_(uplo, &i__1, &z__1, &a[k * a_dim1 + 1], &c__1, &a[
			a_offset], lda, (ftnlen)1);

/*              Store U(k) in column k */

#line 380 "csytf2.f"
		i__1 = k - 1;
#line 380 "csytf2.f"
		cscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 381 "csytf2.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */

#line 395 "csytf2.f"
		if (k > 2) {

#line 397 "csytf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 397 "csytf2.f"
		    d12.r = a[i__1].r, d12.i = a[i__1].i;
#line 398 "csytf2.f"
		    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &d12);
#line 398 "csytf2.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 399 "csytf2.f"
		    z_div(&z__1, &a[k + k * a_dim1], &d12);
#line 399 "csytf2.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 400 "csytf2.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 400 "csytf2.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 400 "csytf2.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 400 "csytf2.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 401 "csytf2.f"
		    z_div(&z__1, &t, &d12);
#line 401 "csytf2.f"
		    d12.r = z__1.r, d12.i = z__1.i;

#line 403 "csytf2.f"
		    for (j = k - 2; j >= 1; --j) {
#line 404 "csytf2.f"
			i__1 = j + (k - 1) * a_dim1;
#line 404 "csytf2.f"
			z__3.r = d11.r * a[i__1].r - d11.i * a[i__1].i, 
				z__3.i = d11.r * a[i__1].i + d11.i * a[i__1]
				.r;
#line 404 "csytf2.f"
			i__2 = j + k * a_dim1;
#line 404 "csytf2.f"
			z__2.r = z__3.r - a[i__2].r, z__2.i = z__3.i - a[i__2]
				.i;
#line 404 "csytf2.f"
			z__1.r = d12.r * z__2.r - d12.i * z__2.i, z__1.i = 
				d12.r * z__2.i + d12.i * z__2.r;
#line 404 "csytf2.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 405 "csytf2.f"
			i__1 = j + k * a_dim1;
#line 405 "csytf2.f"
			z__3.r = d22.r * a[i__1].r - d22.i * a[i__1].i, 
				z__3.i = d22.r * a[i__1].i + d22.i * a[i__1]
				.r;
#line 405 "csytf2.f"
			i__2 = j + (k - 1) * a_dim1;
#line 405 "csytf2.f"
			z__2.r = z__3.r - a[i__2].r, z__2.i = z__3.i - a[i__2]
				.i;
#line 405 "csytf2.f"
			z__1.r = d12.r * z__2.r - d12.i * z__2.i, z__1.i = 
				d12.r * z__2.i + d12.i * z__2.r;
#line 405 "csytf2.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 406 "csytf2.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 407 "csytf2.f"
			    i__1 = i__ + j * a_dim1;
#line 407 "csytf2.f"
			    i__2 = i__ + j * a_dim1;
#line 407 "csytf2.f"
			    i__3 = i__ + k * a_dim1;
#line 407 "csytf2.f"
			    z__3.r = a[i__3].r * wk.r - a[i__3].i * wk.i, 
				    z__3.i = a[i__3].r * wk.i + a[i__3].i * 
				    wk.r;
#line 407 "csytf2.f"
			    z__2.r = a[i__2].r - z__3.r, z__2.i = a[i__2].i - 
				    z__3.i;
#line 407 "csytf2.f"
			    i__4 = i__ + (k - 1) * a_dim1;
#line 407 "csytf2.f"
			    z__4.r = a[i__4].r * wkm1.r - a[i__4].i * wkm1.i, 
				    z__4.i = a[i__4].r * wkm1.i + a[i__4].i * 
				    wkm1.r;
#line 407 "csytf2.f"
			    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - 
				    z__4.i;
#line 407 "csytf2.f"
			    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 409 "csytf2.f"
/* L20: */
#line 409 "csytf2.f"
			}
#line 410 "csytf2.f"
			i__1 = j + k * a_dim1;
#line 410 "csytf2.f"
			a[i__1].r = wk.r, a[i__1].i = wk.i;
#line 411 "csytf2.f"
			i__1 = j + (k - 1) * a_dim1;
#line 411 "csytf2.f"
			a[i__1].r = wkm1.r, a[i__1].i = wkm1.i;
#line 412 "csytf2.f"
/* L30: */
#line 412 "csytf2.f"
		    }

#line 414 "csytf2.f"
		}

#line 416 "csytf2.f"
	    }
#line 417 "csytf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 421 "csytf2.f"
	if (kstep == 1) {
#line 422 "csytf2.f"
	    ipiv[k] = kp;
#line 423 "csytf2.f"
	} else {
#line 424 "csytf2.f"
	    ipiv[k] = -kp;
#line 425 "csytf2.f"
	    ipiv[k - 1] = -kp;
#line 426 "csytf2.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 430 "csytf2.f"
	k -= kstep;
#line 431 "csytf2.f"
	goto L10;

#line 433 "csytf2.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 440 "csytf2.f"
	k = 1;
#line 441 "csytf2.f"
L40:

/*        If K > N, exit from loop */

#line 445 "csytf2.f"
	if (k > *n) {
#line 445 "csytf2.f"
	    goto L70;
#line 445 "csytf2.f"
	}
#line 447 "csytf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 452 "csytf2.f"
	i__1 = k + k * a_dim1;
#line 452 "csytf2.f"
	absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * 
		a_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 458 "csytf2.f"
	if (k < *n) {
#line 459 "csytf2.f"
	    i__1 = *n - k;
#line 459 "csytf2.f"
	    imax = k + icamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 460 "csytf2.f"
	    i__1 = imax + k * a_dim1;
#line 460 "csytf2.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 461 "csytf2.f"
	} else {
#line 462 "csytf2.f"
	    colmax = 0.;
#line 463 "csytf2.f"
	}

#line 465 "csytf2.f"
	if (max(absakk,colmax) == 0. || sisnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 470 "csytf2.f"
	    if (*info == 0) {
#line 470 "csytf2.f"
		*info = k;
#line 470 "csytf2.f"
	    }
#line 472 "csytf2.f"
	    kp = k;
#line 473 "csytf2.f"
	} else {
#line 474 "csytf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 478 "csytf2.f"
		kp = k;
#line 479 "csytf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 484 "csytf2.f"
		i__1 = imax - k;
#line 484 "csytf2.f"
		jmax = k - 1 + icamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 485 "csytf2.f"
		i__1 = imax + jmax * a_dim1;
#line 485 "csytf2.f"
		rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			imax + jmax * a_dim1]), abs(d__2));
#line 486 "csytf2.f"
		if (imax < *n) {
#line 487 "csytf2.f"
		    i__1 = *n - imax;
#line 487 "csytf2.f"
		    jmax = imax + icamax_(&i__1, &a[imax + 1 + imax * a_dim1],
			     &c__1);
/* Computing MAX */
#line 488 "csytf2.f"
		    i__1 = jmax + imax * a_dim1;
#line 488 "csytf2.f"
		    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&a[jmax + imax * a_dim1]), abs(d__2)
			    );
#line 488 "csytf2.f"
		    rowmax = max(d__3,d__4);
#line 489 "csytf2.f"
		}

#line 491 "csytf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 495 "csytf2.f"
		    kp = k;
#line 496 "csytf2.f"
		} else /* if(complicated condition) */ {
#line 496 "csytf2.f"
		    i__1 = imax + imax * a_dim1;
#line 496 "csytf2.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			    imax + imax * a_dim1]), abs(d__2)) >= alpha * 
			    rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 501 "csytf2.f"
			kp = imax;
#line 502 "csytf2.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 507 "csytf2.f"
			kp = imax;
#line 508 "csytf2.f"
			kstep = 2;
#line 509 "csytf2.f"
		    }
#line 509 "csytf2.f"
		}
#line 510 "csytf2.f"
	    }

#line 512 "csytf2.f"
	    kk = k + kstep - 1;
#line 513 "csytf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 518 "csytf2.f"
		if (kp < *n) {
#line 518 "csytf2.f"
		    i__1 = *n - kp;
#line 518 "csytf2.f"
		    cswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 518 "csytf2.f"
		}
#line 520 "csytf2.f"
		i__1 = kp - kk - 1;
#line 520 "csytf2.f"
		cswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 522 "csytf2.f"
		i__1 = kk + kk * a_dim1;
#line 522 "csytf2.f"
		t.r = a[i__1].r, t.i = a[i__1].i;
#line 523 "csytf2.f"
		i__1 = kk + kk * a_dim1;
#line 523 "csytf2.f"
		i__2 = kp + kp * a_dim1;
#line 523 "csytf2.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 524 "csytf2.f"
		i__1 = kp + kp * a_dim1;
#line 524 "csytf2.f"
		a[i__1].r = t.r, a[i__1].i = t.i;
#line 525 "csytf2.f"
		if (kstep == 2) {
#line 526 "csytf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 526 "csytf2.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 527 "csytf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 527 "csytf2.f"
		    i__2 = kp + k * a_dim1;
#line 527 "csytf2.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 528 "csytf2.f"
		    i__1 = kp + k * a_dim1;
#line 528 "csytf2.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 529 "csytf2.f"
		}
#line 530 "csytf2.f"
	    }

/*           Update the trailing submatrix */

#line 534 "csytf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 542 "csytf2.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */

#line 548 "csytf2.f"
		    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 548 "csytf2.f"
		    r1.r = z__1.r, r1.i = z__1.i;
#line 549 "csytf2.f"
		    i__1 = *n - k;
#line 549 "csytf2.f"
		    z__1.r = -r1.r, z__1.i = -r1.i;
#line 549 "csytf2.f"
		    csyr_(uplo, &i__1, &z__1, &a[k + 1 + k * a_dim1], &c__1, &
			    a[k + 1 + (k + 1) * a_dim1], lda, (ftnlen)1);

/*                 Store L(k) in column K */

#line 554 "csytf2.f"
		    i__1 = *n - k;
#line 554 "csytf2.f"
		    cscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 555 "csytf2.f"
		}
#line 556 "csytf2.f"
	    } else {

/*              2-by-2 pivot block D(k) */

#line 560 "csytf2.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T */
/*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 570 "csytf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 570 "csytf2.f"
		    d21.r = a[i__1].r, d21.i = a[i__1].i;
#line 571 "csytf2.f"
		    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &d21);
#line 571 "csytf2.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 572 "csytf2.f"
		    z_div(&z__1, &a[k + k * a_dim1], &d21);
#line 572 "csytf2.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 573 "csytf2.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 573 "csytf2.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 573 "csytf2.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 573 "csytf2.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 574 "csytf2.f"
		    z_div(&z__1, &t, &d21);
#line 574 "csytf2.f"
		    d21.r = z__1.r, d21.i = z__1.i;

#line 576 "csytf2.f"
		    i__1 = *n;
#line 576 "csytf2.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 577 "csytf2.f"
			i__2 = j + k * a_dim1;
#line 577 "csytf2.f"
			z__3.r = d11.r * a[i__2].r - d11.i * a[i__2].i, 
				z__3.i = d11.r * a[i__2].i + d11.i * a[i__2]
				.r;
#line 577 "csytf2.f"
			i__3 = j + (k + 1) * a_dim1;
#line 577 "csytf2.f"
			z__2.r = z__3.r - a[i__3].r, z__2.i = z__3.i - a[i__3]
				.i;
#line 577 "csytf2.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 577 "csytf2.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 578 "csytf2.f"
			i__2 = j + (k + 1) * a_dim1;
#line 578 "csytf2.f"
			z__3.r = d22.r * a[i__2].r - d22.i * a[i__2].i, 
				z__3.i = d22.r * a[i__2].i + d22.i * a[i__2]
				.r;
#line 578 "csytf2.f"
			i__3 = j + k * a_dim1;
#line 578 "csytf2.f"
			z__2.r = z__3.r - a[i__3].r, z__2.i = z__3.i - a[i__3]
				.i;
#line 578 "csytf2.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 578 "csytf2.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;
#line 579 "csytf2.f"
			i__2 = *n;
#line 579 "csytf2.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 580 "csytf2.f"
			    i__3 = i__ + j * a_dim1;
#line 580 "csytf2.f"
			    i__4 = i__ + j * a_dim1;
#line 580 "csytf2.f"
			    i__5 = i__ + k * a_dim1;
#line 580 "csytf2.f"
			    z__3.r = a[i__5].r * wk.r - a[i__5].i * wk.i, 
				    z__3.i = a[i__5].r * wk.i + a[i__5].i * 
				    wk.r;
#line 580 "csytf2.f"
			    z__2.r = a[i__4].r - z__3.r, z__2.i = a[i__4].i - 
				    z__3.i;
#line 580 "csytf2.f"
			    i__6 = i__ + (k + 1) * a_dim1;
#line 580 "csytf2.f"
			    z__4.r = a[i__6].r * wkp1.r - a[i__6].i * wkp1.i, 
				    z__4.i = a[i__6].r * wkp1.i + a[i__6].i * 
				    wkp1.r;
#line 580 "csytf2.f"
			    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - 
				    z__4.i;
#line 580 "csytf2.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 582 "csytf2.f"
/* L50: */
#line 582 "csytf2.f"
			}
#line 583 "csytf2.f"
			i__2 = j + k * a_dim1;
#line 583 "csytf2.f"
			a[i__2].r = wk.r, a[i__2].i = wk.i;
#line 584 "csytf2.f"
			i__2 = j + (k + 1) * a_dim1;
#line 584 "csytf2.f"
			a[i__2].r = wkp1.r, a[i__2].i = wkp1.i;
#line 585 "csytf2.f"
/* L60: */
#line 585 "csytf2.f"
		    }
#line 586 "csytf2.f"
		}
#line 587 "csytf2.f"
	    }
#line 588 "csytf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 592 "csytf2.f"
	if (kstep == 1) {
#line 593 "csytf2.f"
	    ipiv[k] = kp;
#line 594 "csytf2.f"
	} else {
#line 595 "csytf2.f"
	    ipiv[k] = -kp;
#line 596 "csytf2.f"
	    ipiv[k + 1] = -kp;
#line 597 "csytf2.f"
	}

/*        Increase K and return to the start of the main loop */

#line 601 "csytf2.f"
	k += kstep;
#line 602 "csytf2.f"
	goto L40;

#line 604 "csytf2.f"
    }

#line 606 "csytf2.f"
L70:
#line 607 "csytf2.f"
    return 0;

/*     End of CSYTF2 */

} /* csytf2_ */

