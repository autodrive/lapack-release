#line 1 "zhetf2.f"
/* zhetf2.f -- translated by f2c (version 20100827).
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

#line 1 "zhetf2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZHETF2 computes the factorization of a complex Hermitian matrix, using the diagonal pivoting me
thod (unblocked algorithm, calling Level 2 BLAS). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETF2( UPLO, N, A, LDA, IPIV, INFO ) */

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
/* > ZHETF2 computes the factorization of a complex Hermitian matrix A */
/* > using the Bunch-Kaufman diagonal pivoting method: */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \ingroup complex16HEcomputational */

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
/* >  09-29-06 - patch from */
/* >    Bobby Cheng, MathWorks */
/* > */
/* >    Replace l.210 and l.393 */
/* >         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN */
/* >    by */
/* >         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN */
/* > */
/* >  01-01-96 - Based on modifications by */
/* >    J. Lewis, Boeing Computer Services Company */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zhetf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
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
    static integer kk, kp;
    static doublecomplex wk;
    static doublereal tt;
    static doublecomplex wkm1, wkp1;
    static integer imax, jmax;
    extern /* Subroutine */ int zher_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal absakk;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static doublereal colmax;
    extern integer izamax_(integer *, doublecomplex *, integer *);
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

#line 245 "zhetf2.f"
    /* Parameter adjustments */
#line 245 "zhetf2.f"
    a_dim1 = *lda;
#line 245 "zhetf2.f"
    a_offset = 1 + a_dim1;
#line 245 "zhetf2.f"
    a -= a_offset;
#line 245 "zhetf2.f"
    --ipiv;
#line 245 "zhetf2.f"

#line 245 "zhetf2.f"
    /* Function Body */
#line 245 "zhetf2.f"
    *info = 0;
#line 246 "zhetf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 247 "zhetf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 248 "zhetf2.f"
	*info = -1;
#line 249 "zhetf2.f"
    } else if (*n < 0) {
#line 250 "zhetf2.f"
	*info = -2;
#line 251 "zhetf2.f"
    } else if (*lda < max(1,*n)) {
#line 252 "zhetf2.f"
	*info = -4;
#line 253 "zhetf2.f"
    }
#line 254 "zhetf2.f"
    if (*info != 0) {
#line 255 "zhetf2.f"
	i__1 = -(*info);
#line 255 "zhetf2.f"
	xerbla_("ZHETF2", &i__1, (ftnlen)6);
#line 256 "zhetf2.f"
	return 0;
#line 257 "zhetf2.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 261 "zhetf2.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 263 "zhetf2.f"
    if (upper) {

/*        Factorize A as U*D*U**H using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 270 "zhetf2.f"
	k = *n;
#line 271 "zhetf2.f"
L10:

/*        If K < 1, exit from loop */

#line 275 "zhetf2.f"
	if (k < 1) {
#line 275 "zhetf2.f"
	    goto L90;
#line 275 "zhetf2.f"
	}
#line 277 "zhetf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 282 "zhetf2.f"
	i__1 = k + k * a_dim1;
#line 282 "zhetf2.f"
	absakk = (d__1 = a[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 288 "zhetf2.f"
	if (k > 1) {
#line 289 "zhetf2.f"
	    i__1 = k - 1;
#line 289 "zhetf2.f"
	    imax = izamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 290 "zhetf2.f"
	    i__1 = imax + k * a_dim1;
#line 290 "zhetf2.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 291 "zhetf2.f"
	} else {
#line 292 "zhetf2.f"
	    colmax = 0.;
#line 293 "zhetf2.f"
	}

#line 295 "zhetf2.f"
	if (max(absakk,colmax) == 0. || disnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 300 "zhetf2.f"
	    if (*info == 0) {
#line 300 "zhetf2.f"
		*info = k;
#line 300 "zhetf2.f"
	    }
#line 302 "zhetf2.f"
	    kp = k;
#line 303 "zhetf2.f"
	    i__1 = k + k * a_dim1;
#line 303 "zhetf2.f"
	    i__2 = k + k * a_dim1;
#line 303 "zhetf2.f"
	    d__1 = a[i__2].r;
#line 303 "zhetf2.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 304 "zhetf2.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

#line 310 "zhetf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 314 "zhetf2.f"
		kp = k;
#line 315 "zhetf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value. */
/*              Determine only ROWMAX. */

#line 321 "zhetf2.f"
		i__1 = k - imax;
#line 321 "zhetf2.f"
		jmax = imax + izamax_(&i__1, &a[imax + (imax + 1) * a_dim1], 
			lda);
#line 322 "zhetf2.f"
		i__1 = imax + jmax * a_dim1;
#line 322 "zhetf2.f"
		rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			imax + jmax * a_dim1]), abs(d__2));
#line 323 "zhetf2.f"
		if (imax > 1) {
#line 324 "zhetf2.f"
		    i__1 = imax - 1;
#line 324 "zhetf2.f"
		    jmax = izamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
#line 325 "zhetf2.f"
		    i__1 = jmax + imax * a_dim1;
#line 325 "zhetf2.f"
		    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&a[jmax + imax * a_dim1]), abs(d__2)
			    );
#line 325 "zhetf2.f"
		    rowmax = max(d__3,d__4);
#line 326 "zhetf2.f"
		}

#line 328 "zhetf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 332 "zhetf2.f"
		    kp = k;

#line 334 "zhetf2.f"
		} else /* if(complicated condition) */ {
#line 334 "zhetf2.f"
		    i__1 = imax + imax * a_dim1;
#line 334 "zhetf2.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 340 "zhetf2.f"
			kp = imax;
#line 341 "zhetf2.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 346 "zhetf2.f"
			kp = imax;
#line 347 "zhetf2.f"
			kstep = 2;
#line 348 "zhetf2.f"
		    }
#line 348 "zhetf2.f"
		}

#line 350 "zhetf2.f"
	    }

/*           ============================================================ */

#line 354 "zhetf2.f"
	    kk = k - kstep + 1;
#line 355 "zhetf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 360 "zhetf2.f"
		i__1 = kp - 1;
#line 360 "zhetf2.f"
		zswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1],
			 &c__1);
#line 361 "zhetf2.f"
		i__1 = kk - 1;
#line 361 "zhetf2.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 362 "zhetf2.f"
		    d_cnjg(&z__1, &a[j + kk * a_dim1]);
#line 362 "zhetf2.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 363 "zhetf2.f"
		    i__2 = j + kk * a_dim1;
#line 363 "zhetf2.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 363 "zhetf2.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 364 "zhetf2.f"
		    i__2 = kp + j * a_dim1;
#line 364 "zhetf2.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 365 "zhetf2.f"
/* L20: */
#line 365 "zhetf2.f"
		}
#line 366 "zhetf2.f"
		i__1 = kp + kk * a_dim1;
#line 366 "zhetf2.f"
		d_cnjg(&z__1, &a[kp + kk * a_dim1]);
#line 366 "zhetf2.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 367 "zhetf2.f"
		i__1 = kk + kk * a_dim1;
#line 367 "zhetf2.f"
		r1 = a[i__1].r;
#line 368 "zhetf2.f"
		i__1 = kk + kk * a_dim1;
#line 368 "zhetf2.f"
		i__2 = kp + kp * a_dim1;
#line 368 "zhetf2.f"
		d__1 = a[i__2].r;
#line 368 "zhetf2.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 369 "zhetf2.f"
		i__1 = kp + kp * a_dim1;
#line 369 "zhetf2.f"
		a[i__1].r = r1, a[i__1].i = 0.;
#line 370 "zhetf2.f"
		if (kstep == 2) {
#line 371 "zhetf2.f"
		    i__1 = k + k * a_dim1;
#line 371 "zhetf2.f"
		    i__2 = k + k * a_dim1;
#line 371 "zhetf2.f"
		    d__1 = a[i__2].r;
#line 371 "zhetf2.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 372 "zhetf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 372 "zhetf2.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 373 "zhetf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 373 "zhetf2.f"
		    i__2 = kp + k * a_dim1;
#line 373 "zhetf2.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 374 "zhetf2.f"
		    i__1 = kp + k * a_dim1;
#line 374 "zhetf2.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 375 "zhetf2.f"
		}
#line 376 "zhetf2.f"
	    } else {
#line 377 "zhetf2.f"
		i__1 = k + k * a_dim1;
#line 377 "zhetf2.f"
		i__2 = k + k * a_dim1;
#line 377 "zhetf2.f"
		d__1 = a[i__2].r;
#line 377 "zhetf2.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 378 "zhetf2.f"
		if (kstep == 2) {
#line 378 "zhetf2.f"
		    i__1 = k - 1 + (k - 1) * a_dim1;
#line 378 "zhetf2.f"
		    i__2 = k - 1 + (k - 1) * a_dim1;
#line 378 "zhetf2.f"
		    d__1 = a[i__2].r;
#line 378 "zhetf2.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 378 "zhetf2.f"
		}
#line 380 "zhetf2.f"
	    }

/*           Update the leading submatrix */

#line 384 "zhetf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H */

#line 396 "zhetf2.f"
		i__1 = k + k * a_dim1;
#line 396 "zhetf2.f"
		r1 = 1. / a[i__1].r;
#line 397 "zhetf2.f"
		i__1 = k - 1;
#line 397 "zhetf2.f"
		d__1 = -r1;
#line 397 "zhetf2.f"
		zher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &a[
			a_offset], lda, (ftnlen)1);

/*              Store U(k) in column k */

#line 401 "zhetf2.f"
		i__1 = k - 1;
#line 401 "zhetf2.f"
		zdscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 402 "zhetf2.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H */

#line 416 "zhetf2.f"
		if (k > 2) {

#line 418 "zhetf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 418 "zhetf2.f"
		    d__1 = a[i__1].r;
#line 418 "zhetf2.f"
		    d__2 = d_imag(&a[k - 1 + k * a_dim1]);
#line 418 "zhetf2.f"
		    d__ = dlapy2_(&d__1, &d__2);
#line 420 "zhetf2.f"
		    i__1 = k - 1 + (k - 1) * a_dim1;
#line 420 "zhetf2.f"
		    d22 = a[i__1].r / d__;
#line 421 "zhetf2.f"
		    i__1 = k + k * a_dim1;
#line 421 "zhetf2.f"
		    d11 = a[i__1].r / d__;
#line 422 "zhetf2.f"
		    tt = 1. / (d11 * d22 - 1.);
#line 423 "zhetf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 423 "zhetf2.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 423 "zhetf2.f"
		    d12.r = z__1.r, d12.i = z__1.i;
#line 424 "zhetf2.f"
		    d__ = tt / d__;

#line 426 "zhetf2.f"
		    for (j = k - 2; j >= 1; --j) {
#line 427 "zhetf2.f"
			i__1 = j + (k - 1) * a_dim1;
#line 427 "zhetf2.f"
			z__3.r = d11 * a[i__1].r, z__3.i = d11 * a[i__1].i;
#line 427 "zhetf2.f"
			d_cnjg(&z__5, &d12);
#line 427 "zhetf2.f"
			i__2 = j + k * a_dim1;
#line 427 "zhetf2.f"
			z__4.r = z__5.r * a[i__2].r - z__5.i * a[i__2].i, 
				z__4.i = z__5.r * a[i__2].i + z__5.i * a[i__2]
				.r;
#line 427 "zhetf2.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 427 "zhetf2.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 427 "zhetf2.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 429 "zhetf2.f"
			i__1 = j + k * a_dim1;
#line 429 "zhetf2.f"
			z__3.r = d22 * a[i__1].r, z__3.i = d22 * a[i__1].i;
#line 429 "zhetf2.f"
			i__2 = j + (k - 1) * a_dim1;
#line 429 "zhetf2.f"
			z__4.r = d12.r * a[i__2].r - d12.i * a[i__2].i, 
				z__4.i = d12.r * a[i__2].i + d12.i * a[i__2]
				.r;
#line 429 "zhetf2.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 429 "zhetf2.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 429 "zhetf2.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 430 "zhetf2.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 431 "zhetf2.f"
			    i__1 = i__ + j * a_dim1;
#line 431 "zhetf2.f"
			    i__2 = i__ + j * a_dim1;
#line 431 "zhetf2.f"
			    i__3 = i__ + k * a_dim1;
#line 431 "zhetf2.f"
			    d_cnjg(&z__4, &wk);
#line 431 "zhetf2.f"
			    z__3.r = a[i__3].r * z__4.r - a[i__3].i * z__4.i, 
				    z__3.i = a[i__3].r * z__4.i + a[i__3].i * 
				    z__4.r;
#line 431 "zhetf2.f"
			    z__2.r = a[i__2].r - z__3.r, z__2.i = a[i__2].i - 
				    z__3.i;
#line 431 "zhetf2.f"
			    i__4 = i__ + (k - 1) * a_dim1;
#line 431 "zhetf2.f"
			    d_cnjg(&z__6, &wkm1);
#line 431 "zhetf2.f"
			    z__5.r = a[i__4].r * z__6.r - a[i__4].i * z__6.i, 
				    z__5.i = a[i__4].r * z__6.i + a[i__4].i * 
				    z__6.r;
#line 431 "zhetf2.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 431 "zhetf2.f"
			    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 433 "zhetf2.f"
/* L30: */
#line 433 "zhetf2.f"
			}
#line 434 "zhetf2.f"
			i__1 = j + k * a_dim1;
#line 434 "zhetf2.f"
			a[i__1].r = wk.r, a[i__1].i = wk.i;
#line 435 "zhetf2.f"
			i__1 = j + (k - 1) * a_dim1;
#line 435 "zhetf2.f"
			a[i__1].r = wkm1.r, a[i__1].i = wkm1.i;
#line 436 "zhetf2.f"
			i__1 = j + j * a_dim1;
#line 436 "zhetf2.f"
			i__2 = j + j * a_dim1;
#line 436 "zhetf2.f"
			d__1 = a[i__2].r;
#line 436 "zhetf2.f"
			z__1.r = d__1, z__1.i = 0.;
#line 436 "zhetf2.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 437 "zhetf2.f"
/* L40: */
#line 437 "zhetf2.f"
		    }

#line 439 "zhetf2.f"
		}

#line 441 "zhetf2.f"
	    }
#line 442 "zhetf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 446 "zhetf2.f"
	if (kstep == 1) {
#line 447 "zhetf2.f"
	    ipiv[k] = kp;
#line 448 "zhetf2.f"
	} else {
#line 449 "zhetf2.f"
	    ipiv[k] = -kp;
#line 450 "zhetf2.f"
	    ipiv[k - 1] = -kp;
#line 451 "zhetf2.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 455 "zhetf2.f"
	k -= kstep;
#line 456 "zhetf2.f"
	goto L10;

#line 458 "zhetf2.f"
    } else {

/*        Factorize A as L*D*L**H using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 465 "zhetf2.f"
	k = 1;
#line 466 "zhetf2.f"
L50:

/*        If K > N, exit from loop */

#line 470 "zhetf2.f"
	if (k > *n) {
#line 470 "zhetf2.f"
	    goto L90;
#line 470 "zhetf2.f"
	}
#line 472 "zhetf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 477 "zhetf2.f"
	i__1 = k + k * a_dim1;
#line 477 "zhetf2.f"
	absakk = (d__1 = a[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 483 "zhetf2.f"
	if (k < *n) {
#line 484 "zhetf2.f"
	    i__1 = *n - k;
#line 484 "zhetf2.f"
	    imax = k + izamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 485 "zhetf2.f"
	    i__1 = imax + k * a_dim1;
#line 485 "zhetf2.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 486 "zhetf2.f"
	} else {
#line 487 "zhetf2.f"
	    colmax = 0.;
#line 488 "zhetf2.f"
	}

#line 490 "zhetf2.f"
	if (max(absakk,colmax) == 0. || disnan_(&absakk)) {

/*           Column K is zero or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 495 "zhetf2.f"
	    if (*info == 0) {
#line 495 "zhetf2.f"
		*info = k;
#line 495 "zhetf2.f"
	    }
#line 497 "zhetf2.f"
	    kp = k;
#line 498 "zhetf2.f"
	    i__1 = k + k * a_dim1;
#line 498 "zhetf2.f"
	    i__2 = k + k * a_dim1;
#line 498 "zhetf2.f"
	    d__1 = a[i__2].r;
#line 498 "zhetf2.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 499 "zhetf2.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

#line 505 "zhetf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 509 "zhetf2.f"
		kp = k;
#line 510 "zhetf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value. */
/*              Determine only ROWMAX. */

#line 516 "zhetf2.f"
		i__1 = imax - k;
#line 516 "zhetf2.f"
		jmax = k - 1 + izamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 517 "zhetf2.f"
		i__1 = imax + jmax * a_dim1;
#line 517 "zhetf2.f"
		rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			imax + jmax * a_dim1]), abs(d__2));
#line 518 "zhetf2.f"
		if (imax < *n) {
#line 519 "zhetf2.f"
		    i__1 = *n - imax;
#line 519 "zhetf2.f"
		    jmax = imax + izamax_(&i__1, &a[imax + 1 + imax * a_dim1],
			     &c__1);
/* Computing MAX */
#line 520 "zhetf2.f"
		    i__1 = jmax + imax * a_dim1;
#line 520 "zhetf2.f"
		    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&a[jmax + imax * a_dim1]), abs(d__2)
			    );
#line 520 "zhetf2.f"
		    rowmax = max(d__3,d__4);
#line 521 "zhetf2.f"
		}

#line 523 "zhetf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 527 "zhetf2.f"
		    kp = k;

#line 529 "zhetf2.f"
		} else /* if(complicated condition) */ {
#line 529 "zhetf2.f"
		    i__1 = imax + imax * a_dim1;
#line 529 "zhetf2.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 535 "zhetf2.f"
			kp = imax;
#line 536 "zhetf2.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 541 "zhetf2.f"
			kp = imax;
#line 542 "zhetf2.f"
			kstep = 2;
#line 543 "zhetf2.f"
		    }
#line 543 "zhetf2.f"
		}

#line 545 "zhetf2.f"
	    }

/*           ============================================================ */

#line 549 "zhetf2.f"
	    kk = k + kstep - 1;
#line 550 "zhetf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 555 "zhetf2.f"
		if (kp < *n) {
#line 555 "zhetf2.f"
		    i__1 = *n - kp;
#line 555 "zhetf2.f"
		    zswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 555 "zhetf2.f"
		}
#line 557 "zhetf2.f"
		i__1 = kp - 1;
#line 557 "zhetf2.f"
		for (j = kk + 1; j <= i__1; ++j) {
#line 558 "zhetf2.f"
		    d_cnjg(&z__1, &a[j + kk * a_dim1]);
#line 558 "zhetf2.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 559 "zhetf2.f"
		    i__2 = j + kk * a_dim1;
#line 559 "zhetf2.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 559 "zhetf2.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 560 "zhetf2.f"
		    i__2 = kp + j * a_dim1;
#line 560 "zhetf2.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 561 "zhetf2.f"
/* L60: */
#line 561 "zhetf2.f"
		}
#line 562 "zhetf2.f"
		i__1 = kp + kk * a_dim1;
#line 562 "zhetf2.f"
		d_cnjg(&z__1, &a[kp + kk * a_dim1]);
#line 562 "zhetf2.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 563 "zhetf2.f"
		i__1 = kk + kk * a_dim1;
#line 563 "zhetf2.f"
		r1 = a[i__1].r;
#line 564 "zhetf2.f"
		i__1 = kk + kk * a_dim1;
#line 564 "zhetf2.f"
		i__2 = kp + kp * a_dim1;
#line 564 "zhetf2.f"
		d__1 = a[i__2].r;
#line 564 "zhetf2.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 565 "zhetf2.f"
		i__1 = kp + kp * a_dim1;
#line 565 "zhetf2.f"
		a[i__1].r = r1, a[i__1].i = 0.;
#line 566 "zhetf2.f"
		if (kstep == 2) {
#line 567 "zhetf2.f"
		    i__1 = k + k * a_dim1;
#line 567 "zhetf2.f"
		    i__2 = k + k * a_dim1;
#line 567 "zhetf2.f"
		    d__1 = a[i__2].r;
#line 567 "zhetf2.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 568 "zhetf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 568 "zhetf2.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 569 "zhetf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 569 "zhetf2.f"
		    i__2 = kp + k * a_dim1;
#line 569 "zhetf2.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 570 "zhetf2.f"
		    i__1 = kp + k * a_dim1;
#line 570 "zhetf2.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 571 "zhetf2.f"
		}
#line 572 "zhetf2.f"
	    } else {
#line 573 "zhetf2.f"
		i__1 = k + k * a_dim1;
#line 573 "zhetf2.f"
		i__2 = k + k * a_dim1;
#line 573 "zhetf2.f"
		d__1 = a[i__2].r;
#line 573 "zhetf2.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 574 "zhetf2.f"
		if (kstep == 2) {
#line 574 "zhetf2.f"
		    i__1 = k + 1 + (k + 1) * a_dim1;
#line 574 "zhetf2.f"
		    i__2 = k + 1 + (k + 1) * a_dim1;
#line 574 "zhetf2.f"
		    d__1 = a[i__2].r;
#line 574 "zhetf2.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 574 "zhetf2.f"
		}
#line 576 "zhetf2.f"
	    }

/*           Update the trailing submatrix */

#line 580 "zhetf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 588 "zhetf2.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H */

#line 594 "zhetf2.f"
		    i__1 = k + k * a_dim1;
#line 594 "zhetf2.f"
		    r1 = 1. / a[i__1].r;
#line 595 "zhetf2.f"
		    i__1 = *n - k;
#line 595 "zhetf2.f"
		    d__1 = -r1;
#line 595 "zhetf2.f"
		    zher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &c__1, &
			    a[k + 1 + (k + 1) * a_dim1], lda, (ftnlen)1);

/*                 Store L(k) in column K */

#line 600 "zhetf2.f"
		    i__1 = *n - k;
#line 600 "zhetf2.f"
		    zdscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 601 "zhetf2.f"
		}
#line 602 "zhetf2.f"
	    } else {

/*              2-by-2 pivot block D(k) */

#line 606 "zhetf2.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H */
/*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 616 "zhetf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 616 "zhetf2.f"
		    d__1 = a[i__1].r;
#line 616 "zhetf2.f"
		    d__2 = d_imag(&a[k + 1 + k * a_dim1]);
#line 616 "zhetf2.f"
		    d__ = dlapy2_(&d__1, &d__2);
#line 618 "zhetf2.f"
		    i__1 = k + 1 + (k + 1) * a_dim1;
#line 618 "zhetf2.f"
		    d11 = a[i__1].r / d__;
#line 619 "zhetf2.f"
		    i__1 = k + k * a_dim1;
#line 619 "zhetf2.f"
		    d22 = a[i__1].r / d__;
#line 620 "zhetf2.f"
		    tt = 1. / (d11 * d22 - 1.);
#line 621 "zhetf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 621 "zhetf2.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 621 "zhetf2.f"
		    d21.r = z__1.r, d21.i = z__1.i;
#line 622 "zhetf2.f"
		    d__ = tt / d__;

#line 624 "zhetf2.f"
		    i__1 = *n;
#line 624 "zhetf2.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 625 "zhetf2.f"
			i__2 = j + k * a_dim1;
#line 625 "zhetf2.f"
			z__3.r = d11 * a[i__2].r, z__3.i = d11 * a[i__2].i;
#line 625 "zhetf2.f"
			i__3 = j + (k + 1) * a_dim1;
#line 625 "zhetf2.f"
			z__4.r = d21.r * a[i__3].r - d21.i * a[i__3].i, 
				z__4.i = d21.r * a[i__3].i + d21.i * a[i__3]
				.r;
#line 625 "zhetf2.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 625 "zhetf2.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 625 "zhetf2.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 626 "zhetf2.f"
			i__2 = j + (k + 1) * a_dim1;
#line 626 "zhetf2.f"
			z__3.r = d22 * a[i__2].r, z__3.i = d22 * a[i__2].i;
#line 626 "zhetf2.f"
			d_cnjg(&z__5, &d21);
#line 626 "zhetf2.f"
			i__3 = j + k * a_dim1;
#line 626 "zhetf2.f"
			z__4.r = z__5.r * a[i__3].r - z__5.i * a[i__3].i, 
				z__4.i = z__5.r * a[i__3].i + z__5.i * a[i__3]
				.r;
#line 626 "zhetf2.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 626 "zhetf2.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 626 "zhetf2.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;
#line 628 "zhetf2.f"
			i__2 = *n;
#line 628 "zhetf2.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 629 "zhetf2.f"
			    i__3 = i__ + j * a_dim1;
#line 629 "zhetf2.f"
			    i__4 = i__ + j * a_dim1;
#line 629 "zhetf2.f"
			    i__5 = i__ + k * a_dim1;
#line 629 "zhetf2.f"
			    d_cnjg(&z__4, &wk);
#line 629 "zhetf2.f"
			    z__3.r = a[i__5].r * z__4.r - a[i__5].i * z__4.i, 
				    z__3.i = a[i__5].r * z__4.i + a[i__5].i * 
				    z__4.r;
#line 629 "zhetf2.f"
			    z__2.r = a[i__4].r - z__3.r, z__2.i = a[i__4].i - 
				    z__3.i;
#line 629 "zhetf2.f"
			    i__6 = i__ + (k + 1) * a_dim1;
#line 629 "zhetf2.f"
			    d_cnjg(&z__6, &wkp1);
#line 629 "zhetf2.f"
			    z__5.r = a[i__6].r * z__6.r - a[i__6].i * z__6.i, 
				    z__5.i = a[i__6].r * z__6.i + a[i__6].i * 
				    z__6.r;
#line 629 "zhetf2.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 629 "zhetf2.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 631 "zhetf2.f"
/* L70: */
#line 631 "zhetf2.f"
			}
#line 632 "zhetf2.f"
			i__2 = j + k * a_dim1;
#line 632 "zhetf2.f"
			a[i__2].r = wk.r, a[i__2].i = wk.i;
#line 633 "zhetf2.f"
			i__2 = j + (k + 1) * a_dim1;
#line 633 "zhetf2.f"
			a[i__2].r = wkp1.r, a[i__2].i = wkp1.i;
#line 634 "zhetf2.f"
			i__2 = j + j * a_dim1;
#line 634 "zhetf2.f"
			i__3 = j + j * a_dim1;
#line 634 "zhetf2.f"
			d__1 = a[i__3].r;
#line 634 "zhetf2.f"
			z__1.r = d__1, z__1.i = 0.;
#line 634 "zhetf2.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 635 "zhetf2.f"
/* L80: */
#line 635 "zhetf2.f"
		    }
#line 636 "zhetf2.f"
		}
#line 637 "zhetf2.f"
	    }
#line 638 "zhetf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 642 "zhetf2.f"
	if (kstep == 1) {
#line 643 "zhetf2.f"
	    ipiv[k] = kp;
#line 644 "zhetf2.f"
	} else {
#line 645 "zhetf2.f"
	    ipiv[k] = -kp;
#line 646 "zhetf2.f"
	    ipiv[k + 1] = -kp;
#line 647 "zhetf2.f"
	}

/*        Increase K and return to the start of the main loop */

#line 651 "zhetf2.f"
	k += kstep;
#line 652 "zhetf2.f"
	goto L50;

#line 654 "zhetf2.f"
    }

#line 656 "zhetf2.f"
L90:
#line 657 "zhetf2.f"
    return 0;

/*     End of ZHETF2 */

} /* zhetf2_ */

