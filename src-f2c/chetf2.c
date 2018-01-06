#line 1 "chetf2.f"
/* chetf2.f -- translated by f2c (version 20100827).
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

#line 1 "chetf2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CHETF2 computes the factorization of a complex Hermitian matrix, using the diagonal pivoting me
thod (unblocked algorithm calling Level 2 BLAS). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETF2( UPLO, N, A, LDA, IPIV, INFO ) */

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
/* > CHETF2 computes the factorization of a complex Hermitian matrix A */
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

/* > \ingroup complexHEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  09-29-06 - patch from */
/* >    Bobby Cheng, MathWorks */
/* > */
/* >    Replace l.210 and l.392 */
/* >         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN */
/* >    by */
/* >         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN */
/* > */
/* >  01-01-96 - Based on modifications by */
/* >    J. Lewis, Boeing Computer Services Company */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
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
/* > */
/*  ===================================================================== */
/* Subroutine */ int chetf2_(char *uplo, integer *n, doublecomplex *a, 
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
    extern /* Subroutine */ int cher_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static integer imax, jmax;
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer kstep;
    static logical upper;
    extern doublereal slapy2_(doublereal *, doublereal *);
    static doublereal absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 240 "chetf2.f"
    /* Parameter adjustments */
#line 240 "chetf2.f"
    a_dim1 = *lda;
#line 240 "chetf2.f"
    a_offset = 1 + a_dim1;
#line 240 "chetf2.f"
    a -= a_offset;
#line 240 "chetf2.f"
    --ipiv;
#line 240 "chetf2.f"

#line 240 "chetf2.f"
    /* Function Body */
#line 240 "chetf2.f"
    *info = 0;
#line 241 "chetf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 242 "chetf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 243 "chetf2.f"
	*info = -1;
#line 244 "chetf2.f"
    } else if (*n < 0) {
#line 245 "chetf2.f"
	*info = -2;
#line 246 "chetf2.f"
    } else if (*lda < max(1,*n)) {
#line 247 "chetf2.f"
	*info = -4;
#line 248 "chetf2.f"
    }
#line 249 "chetf2.f"
    if (*info != 0) {
#line 250 "chetf2.f"
	i__1 = -(*info);
#line 250 "chetf2.f"
	xerbla_("CHETF2", &i__1, (ftnlen)6);
#line 251 "chetf2.f"
	return 0;
#line 252 "chetf2.f"
    }

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 256 "chetf2.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 258 "chetf2.f"
    if (upper) {

/*        Factorize A as U*D*U**H using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2 */

#line 265 "chetf2.f"
	k = *n;
#line 266 "chetf2.f"
L10:

/*        If K < 1, exit from loop */

#line 270 "chetf2.f"
	if (k < 1) {
#line 270 "chetf2.f"
	    goto L90;
#line 270 "chetf2.f"
	}
#line 272 "chetf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 277 "chetf2.f"
	i__1 = k + k * a_dim1;
#line 277 "chetf2.f"
	absakk = (d__1 = a[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 283 "chetf2.f"
	if (k > 1) {
#line 284 "chetf2.f"
	    i__1 = k - 1;
#line 284 "chetf2.f"
	    imax = icamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
#line 285 "chetf2.f"
	    i__1 = imax + k * a_dim1;
#line 285 "chetf2.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 286 "chetf2.f"
	} else {
#line 287 "chetf2.f"
	    colmax = 0.;
#line 288 "chetf2.f"
	}

#line 290 "chetf2.f"
	if (max(absakk,colmax) == 0. || sisnan_(&absakk)) {

/*           Column K is or underflow, or contains a NaN: */
/*           set INFO and continue */

#line 295 "chetf2.f"
	    if (*info == 0) {
#line 295 "chetf2.f"
		*info = k;
#line 295 "chetf2.f"
	    }
#line 297 "chetf2.f"
	    kp = k;
#line 298 "chetf2.f"
	    i__1 = k + k * a_dim1;
#line 298 "chetf2.f"
	    i__2 = k + k * a_dim1;
#line 298 "chetf2.f"
	    d__1 = a[i__2].r;
#line 298 "chetf2.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 299 "chetf2.f"
	} else {
#line 300 "chetf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 304 "chetf2.f"
		kp = k;
#line 305 "chetf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 310 "chetf2.f"
		i__1 = k - imax;
#line 310 "chetf2.f"
		jmax = imax + icamax_(&i__1, &a[imax + (imax + 1) * a_dim1], 
			lda);
#line 311 "chetf2.f"
		i__1 = imax + jmax * a_dim1;
#line 311 "chetf2.f"
		rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			imax + jmax * a_dim1]), abs(d__2));
#line 312 "chetf2.f"
		if (imax > 1) {
#line 313 "chetf2.f"
		    i__1 = imax - 1;
#line 313 "chetf2.f"
		    jmax = icamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
#line 314 "chetf2.f"
		    i__1 = jmax + imax * a_dim1;
#line 314 "chetf2.f"
		    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&a[jmax + imax * a_dim1]), abs(d__2)
			    );
#line 314 "chetf2.f"
		    rowmax = max(d__3,d__4);
#line 315 "chetf2.f"
		}

#line 317 "chetf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 321 "chetf2.f"
		    kp = k;
#line 322 "chetf2.f"
		} else /* if(complicated condition) */ {
#line 322 "chetf2.f"
		    i__1 = imax + imax * a_dim1;
#line 322 "chetf2.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 328 "chetf2.f"
			kp = imax;
#line 329 "chetf2.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 334 "chetf2.f"
			kp = imax;
#line 335 "chetf2.f"
			kstep = 2;
#line 336 "chetf2.f"
		    }
#line 336 "chetf2.f"
		}
#line 337 "chetf2.f"
	    }

#line 339 "chetf2.f"
	    kk = k - kstep + 1;
#line 340 "chetf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the leading */
/*              submatrix A(1:k,1:k) */

#line 345 "chetf2.f"
		i__1 = kp - 1;
#line 345 "chetf2.f"
		cswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1],
			 &c__1);
#line 346 "chetf2.f"
		i__1 = kk - 1;
#line 346 "chetf2.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 347 "chetf2.f"
		    d_cnjg(&z__1, &a[j + kk * a_dim1]);
#line 347 "chetf2.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 348 "chetf2.f"
		    i__2 = j + kk * a_dim1;
#line 348 "chetf2.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 348 "chetf2.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 349 "chetf2.f"
		    i__2 = kp + j * a_dim1;
#line 349 "chetf2.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 350 "chetf2.f"
/* L20: */
#line 350 "chetf2.f"
		}
#line 351 "chetf2.f"
		i__1 = kp + kk * a_dim1;
#line 351 "chetf2.f"
		d_cnjg(&z__1, &a[kp + kk * a_dim1]);
#line 351 "chetf2.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 352 "chetf2.f"
		i__1 = kk + kk * a_dim1;
#line 352 "chetf2.f"
		r1 = a[i__1].r;
#line 353 "chetf2.f"
		i__1 = kk + kk * a_dim1;
#line 353 "chetf2.f"
		i__2 = kp + kp * a_dim1;
#line 353 "chetf2.f"
		d__1 = a[i__2].r;
#line 353 "chetf2.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 354 "chetf2.f"
		i__1 = kp + kp * a_dim1;
#line 354 "chetf2.f"
		a[i__1].r = r1, a[i__1].i = 0.;
#line 355 "chetf2.f"
		if (kstep == 2) {
#line 356 "chetf2.f"
		    i__1 = k + k * a_dim1;
#line 356 "chetf2.f"
		    i__2 = k + k * a_dim1;
#line 356 "chetf2.f"
		    d__1 = a[i__2].r;
#line 356 "chetf2.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 357 "chetf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 357 "chetf2.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 358 "chetf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 358 "chetf2.f"
		    i__2 = kp + k * a_dim1;
#line 358 "chetf2.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 359 "chetf2.f"
		    i__1 = kp + k * a_dim1;
#line 359 "chetf2.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 360 "chetf2.f"
		}
#line 361 "chetf2.f"
	    } else {
#line 362 "chetf2.f"
		i__1 = k + k * a_dim1;
#line 362 "chetf2.f"
		i__2 = k + k * a_dim1;
#line 362 "chetf2.f"
		d__1 = a[i__2].r;
#line 362 "chetf2.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 363 "chetf2.f"
		if (kstep == 2) {
#line 363 "chetf2.f"
		    i__1 = k - 1 + (k - 1) * a_dim1;
#line 363 "chetf2.f"
		    i__2 = k - 1 + (k - 1) * a_dim1;
#line 363 "chetf2.f"
		    d__1 = a[i__2].r;
#line 363 "chetf2.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 363 "chetf2.f"
		}
#line 365 "chetf2.f"
	    }

/*           Update the leading submatrix */

#line 369 "chetf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Perform a rank-1 update of A(1:k-1,1:k-1) as */

/*              A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H */

#line 381 "chetf2.f"
		i__1 = k + k * a_dim1;
#line 381 "chetf2.f"
		r1 = 1. / a[i__1].r;
#line 382 "chetf2.f"
		i__1 = k - 1;
#line 382 "chetf2.f"
		d__1 = -r1;
#line 382 "chetf2.f"
		cher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &a[
			a_offset], lda, (ftnlen)1);

/*              Store U(k) in column k */

#line 386 "chetf2.f"
		i__1 = k - 1;
#line 386 "chetf2.f"
		csscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 387 "chetf2.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k-1 now hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Perform a rank-2 update of A(1:k-2,1:k-2) as */

/*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H */
/*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H */

#line 401 "chetf2.f"
		if (k > 2) {

#line 403 "chetf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 403 "chetf2.f"
		    d__1 = a[i__1].r;
#line 403 "chetf2.f"
		    d__2 = d_imag(&a[k - 1 + k * a_dim1]);
#line 403 "chetf2.f"
		    d__ = slapy2_(&d__1, &d__2);
#line 405 "chetf2.f"
		    i__1 = k - 1 + (k - 1) * a_dim1;
#line 405 "chetf2.f"
		    d22 = a[i__1].r / d__;
#line 406 "chetf2.f"
		    i__1 = k + k * a_dim1;
#line 406 "chetf2.f"
		    d11 = a[i__1].r / d__;
#line 407 "chetf2.f"
		    tt = 1. / (d11 * d22 - 1.);
#line 408 "chetf2.f"
		    i__1 = k - 1 + k * a_dim1;
#line 408 "chetf2.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 408 "chetf2.f"
		    d12.r = z__1.r, d12.i = z__1.i;
#line 409 "chetf2.f"
		    d__ = tt / d__;

#line 411 "chetf2.f"
		    for (j = k - 2; j >= 1; --j) {
#line 412 "chetf2.f"
			i__1 = j + (k - 1) * a_dim1;
#line 412 "chetf2.f"
			z__3.r = d11 * a[i__1].r, z__3.i = d11 * a[i__1].i;
#line 412 "chetf2.f"
			d_cnjg(&z__5, &d12);
#line 412 "chetf2.f"
			i__2 = j + k * a_dim1;
#line 412 "chetf2.f"
			z__4.r = z__5.r * a[i__2].r - z__5.i * a[i__2].i, 
				z__4.i = z__5.r * a[i__2].i + z__5.i * a[i__2]
				.r;
#line 412 "chetf2.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 412 "chetf2.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 412 "chetf2.f"
			wkm1.r = z__1.r, wkm1.i = z__1.i;
#line 413 "chetf2.f"
			i__1 = j + k * a_dim1;
#line 413 "chetf2.f"
			z__3.r = d22 * a[i__1].r, z__3.i = d22 * a[i__1].i;
#line 413 "chetf2.f"
			i__2 = j + (k - 1) * a_dim1;
#line 413 "chetf2.f"
			z__4.r = d12.r * a[i__2].r - d12.i * a[i__2].i, 
				z__4.i = d12.r * a[i__2].i + d12.i * a[i__2]
				.r;
#line 413 "chetf2.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 413 "chetf2.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 413 "chetf2.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 414 "chetf2.f"
			for (i__ = j; i__ >= 1; --i__) {
#line 415 "chetf2.f"
			    i__1 = i__ + j * a_dim1;
#line 415 "chetf2.f"
			    i__2 = i__ + j * a_dim1;
#line 415 "chetf2.f"
			    i__3 = i__ + k * a_dim1;
#line 415 "chetf2.f"
			    d_cnjg(&z__4, &wk);
#line 415 "chetf2.f"
			    z__3.r = a[i__3].r * z__4.r - a[i__3].i * z__4.i, 
				    z__3.i = a[i__3].r * z__4.i + a[i__3].i * 
				    z__4.r;
#line 415 "chetf2.f"
			    z__2.r = a[i__2].r - z__3.r, z__2.i = a[i__2].i - 
				    z__3.i;
#line 415 "chetf2.f"
			    i__4 = i__ + (k - 1) * a_dim1;
#line 415 "chetf2.f"
			    d_cnjg(&z__6, &wkm1);
#line 415 "chetf2.f"
			    z__5.r = a[i__4].r * z__6.r - a[i__4].i * z__6.i, 
				    z__5.i = a[i__4].r * z__6.i + a[i__4].i * 
				    z__6.r;
#line 415 "chetf2.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 415 "chetf2.f"
			    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 417 "chetf2.f"
/* L30: */
#line 417 "chetf2.f"
			}
#line 418 "chetf2.f"
			i__1 = j + k * a_dim1;
#line 418 "chetf2.f"
			a[i__1].r = wk.r, a[i__1].i = wk.i;
#line 419 "chetf2.f"
			i__1 = j + (k - 1) * a_dim1;
#line 419 "chetf2.f"
			a[i__1].r = wkm1.r, a[i__1].i = wkm1.i;
#line 420 "chetf2.f"
			i__1 = j + j * a_dim1;
#line 420 "chetf2.f"
			i__2 = j + j * a_dim1;
#line 420 "chetf2.f"
			d__1 = a[i__2].r;
#line 420 "chetf2.f"
			z__1.r = d__1, z__1.i = 0.;
#line 420 "chetf2.f"
			a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 421 "chetf2.f"
/* L40: */
#line 421 "chetf2.f"
		    }

#line 423 "chetf2.f"
		}

#line 425 "chetf2.f"
	    }
#line 426 "chetf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 430 "chetf2.f"
	if (kstep == 1) {
#line 431 "chetf2.f"
	    ipiv[k] = kp;
#line 432 "chetf2.f"
	} else {
#line 433 "chetf2.f"
	    ipiv[k] = -kp;
#line 434 "chetf2.f"
	    ipiv[k - 1] = -kp;
#line 435 "chetf2.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 439 "chetf2.f"
	k -= kstep;
#line 440 "chetf2.f"
	goto L10;

#line 442 "chetf2.f"
    } else {

/*        Factorize A as L*D*L**H using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2 */

#line 449 "chetf2.f"
	k = 1;
#line 450 "chetf2.f"
L50:

/*        If K > N, exit from loop */

#line 454 "chetf2.f"
	if (k > *n) {
#line 454 "chetf2.f"
	    goto L90;
#line 454 "chetf2.f"
	}
#line 456 "chetf2.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 461 "chetf2.f"
	i__1 = k + k * a_dim1;
#line 461 "chetf2.f"
	absakk = (d__1 = a[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 467 "chetf2.f"
	if (k < *n) {
#line 468 "chetf2.f"
	    i__1 = *n - k;
#line 468 "chetf2.f"
	    imax = k + icamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
#line 469 "chetf2.f"
	    i__1 = imax + k * a_dim1;
#line 469 "chetf2.f"
	    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + 
		    k * a_dim1]), abs(d__2));
#line 470 "chetf2.f"
	} else {
#line 471 "chetf2.f"
	    colmax = 0.;
#line 472 "chetf2.f"
	}

#line 474 "chetf2.f"
	if (max(absakk,colmax) == 0. || sisnan_(&absakk)) {

/*           Column K is zero or underflow, contains a NaN: */
/*           set INFO and continue */

#line 479 "chetf2.f"
	    if (*info == 0) {
#line 479 "chetf2.f"
		*info = k;
#line 479 "chetf2.f"
	    }
#line 481 "chetf2.f"
	    kp = k;
#line 482 "chetf2.f"
	    i__1 = k + k * a_dim1;
#line 482 "chetf2.f"
	    i__2 = k + k * a_dim1;
#line 482 "chetf2.f"
	    d__1 = a[i__2].r;
#line 482 "chetf2.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 483 "chetf2.f"
	} else {
#line 484 "chetf2.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 488 "chetf2.f"
		kp = k;
#line 489 "chetf2.f"
	    } else {

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 494 "chetf2.f"
		i__1 = imax - k;
#line 494 "chetf2.f"
		jmax = k - 1 + icamax_(&i__1, &a[imax + k * a_dim1], lda);
#line 495 "chetf2.f"
		i__1 = imax + jmax * a_dim1;
#line 495 "chetf2.f"
		rowmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
			imax + jmax * a_dim1]), abs(d__2));
#line 496 "chetf2.f"
		if (imax < *n) {
#line 497 "chetf2.f"
		    i__1 = *n - imax;
#line 497 "chetf2.f"
		    jmax = imax + icamax_(&i__1, &a[imax + 1 + imax * a_dim1],
			     &c__1);
/* Computing MAX */
#line 498 "chetf2.f"
		    i__1 = jmax + imax * a_dim1;
#line 498 "chetf2.f"
		    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&a[jmax + imax * a_dim1]), abs(d__2)
			    );
#line 498 "chetf2.f"
		    rowmax = max(d__3,d__4);
#line 499 "chetf2.f"
		}

#line 501 "chetf2.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 505 "chetf2.f"
		    kp = k;
#line 506 "chetf2.f"
		} else /* if(complicated condition) */ {
#line 506 "chetf2.f"
		    i__1 = imax + imax * a_dim1;
#line 506 "chetf2.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 512 "chetf2.f"
			kp = imax;
#line 513 "chetf2.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 518 "chetf2.f"
			kp = imax;
#line 519 "chetf2.f"
			kstep = 2;
#line 520 "chetf2.f"
		    }
#line 520 "chetf2.f"
		}
#line 521 "chetf2.f"
	    }

#line 523 "chetf2.f"
	    kk = k + kstep - 1;
#line 524 "chetf2.f"
	    if (kp != kk) {

/*              Interchange rows and columns KK and KP in the trailing */
/*              submatrix A(k:n,k:n) */

#line 529 "chetf2.f"
		if (kp < *n) {
#line 529 "chetf2.f"
		    i__1 = *n - kp;
#line 529 "chetf2.f"
		    cswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 529 "chetf2.f"
		}
#line 531 "chetf2.f"
		i__1 = kp - 1;
#line 531 "chetf2.f"
		for (j = kk + 1; j <= i__1; ++j) {
#line 532 "chetf2.f"
		    d_cnjg(&z__1, &a[j + kk * a_dim1]);
#line 532 "chetf2.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 533 "chetf2.f"
		    i__2 = j + kk * a_dim1;
#line 533 "chetf2.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 533 "chetf2.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 534 "chetf2.f"
		    i__2 = kp + j * a_dim1;
#line 534 "chetf2.f"
		    a[i__2].r = t.r, a[i__2].i = t.i;
#line 535 "chetf2.f"
/* L60: */
#line 535 "chetf2.f"
		}
#line 536 "chetf2.f"
		i__1 = kp + kk * a_dim1;
#line 536 "chetf2.f"
		d_cnjg(&z__1, &a[kp + kk * a_dim1]);
#line 536 "chetf2.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 537 "chetf2.f"
		i__1 = kk + kk * a_dim1;
#line 537 "chetf2.f"
		r1 = a[i__1].r;
#line 538 "chetf2.f"
		i__1 = kk + kk * a_dim1;
#line 538 "chetf2.f"
		i__2 = kp + kp * a_dim1;
#line 538 "chetf2.f"
		d__1 = a[i__2].r;
#line 538 "chetf2.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 539 "chetf2.f"
		i__1 = kp + kp * a_dim1;
#line 539 "chetf2.f"
		a[i__1].r = r1, a[i__1].i = 0.;
#line 540 "chetf2.f"
		if (kstep == 2) {
#line 541 "chetf2.f"
		    i__1 = k + k * a_dim1;
#line 541 "chetf2.f"
		    i__2 = k + k * a_dim1;
#line 541 "chetf2.f"
		    d__1 = a[i__2].r;
#line 541 "chetf2.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 542 "chetf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 542 "chetf2.f"
		    t.r = a[i__1].r, t.i = a[i__1].i;
#line 543 "chetf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 543 "chetf2.f"
		    i__2 = kp + k * a_dim1;
#line 543 "chetf2.f"
		    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 544 "chetf2.f"
		    i__1 = kp + k * a_dim1;
#line 544 "chetf2.f"
		    a[i__1].r = t.r, a[i__1].i = t.i;
#line 545 "chetf2.f"
		}
#line 546 "chetf2.f"
	    } else {
#line 547 "chetf2.f"
		i__1 = k + k * a_dim1;
#line 547 "chetf2.f"
		i__2 = k + k * a_dim1;
#line 547 "chetf2.f"
		d__1 = a[i__2].r;
#line 547 "chetf2.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 548 "chetf2.f"
		if (kstep == 2) {
#line 548 "chetf2.f"
		    i__1 = k + 1 + (k + 1) * a_dim1;
#line 548 "chetf2.f"
		    i__2 = k + 1 + (k + 1) * a_dim1;
#line 548 "chetf2.f"
		    d__1 = a[i__2].r;
#line 548 "chetf2.f"
		    a[i__1].r = d__1, a[i__1].i = 0.;
#line 548 "chetf2.f"
		}
#line 550 "chetf2.f"
	    }

/*           Update the trailing submatrix */

#line 554 "chetf2.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

#line 562 "chetf2.f"
		if (k < *n) {

/*                 Perform a rank-1 update of A(k+1:n,k+1:n) as */

/*                 A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H */

#line 568 "chetf2.f"
		    i__1 = k + k * a_dim1;
#line 568 "chetf2.f"
		    r1 = 1. / a[i__1].r;
#line 569 "chetf2.f"
		    i__1 = *n - k;
#line 569 "chetf2.f"
		    d__1 = -r1;
#line 569 "chetf2.f"
		    cher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &c__1, &
			    a[k + 1 + (k + 1) * a_dim1], lda, (ftnlen)1);

/*                 Store L(k) in column K */

#line 574 "chetf2.f"
		    i__1 = *n - k;
#line 574 "chetf2.f"
		    csscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 575 "chetf2.f"
		}
#line 576 "chetf2.f"
	    } else {

/*              2-by-2 pivot block D(k) */

#line 580 "chetf2.f"
		if (k < *n - 1) {

/*                 Perform a rank-2 update of A(k+2:n,k+2:n) as */

/*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H */
/*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H */

/*                 where L(k) and L(k+1) are the k-th and (k+1)-th */
/*                 columns of L */

#line 590 "chetf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 590 "chetf2.f"
		    d__1 = a[i__1].r;
#line 590 "chetf2.f"
		    d__2 = d_imag(&a[k + 1 + k * a_dim1]);
#line 590 "chetf2.f"
		    d__ = slapy2_(&d__1, &d__2);
#line 592 "chetf2.f"
		    i__1 = k + 1 + (k + 1) * a_dim1;
#line 592 "chetf2.f"
		    d11 = a[i__1].r / d__;
#line 593 "chetf2.f"
		    i__1 = k + k * a_dim1;
#line 593 "chetf2.f"
		    d22 = a[i__1].r / d__;
#line 594 "chetf2.f"
		    tt = 1. / (d11 * d22 - 1.);
#line 595 "chetf2.f"
		    i__1 = k + 1 + k * a_dim1;
#line 595 "chetf2.f"
		    z__1.r = a[i__1].r / d__, z__1.i = a[i__1].i / d__;
#line 595 "chetf2.f"
		    d21.r = z__1.r, d21.i = z__1.i;
#line 596 "chetf2.f"
		    d__ = tt / d__;

#line 598 "chetf2.f"
		    i__1 = *n;
#line 598 "chetf2.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 599 "chetf2.f"
			i__2 = j + k * a_dim1;
#line 599 "chetf2.f"
			z__3.r = d11 * a[i__2].r, z__3.i = d11 * a[i__2].i;
#line 599 "chetf2.f"
			i__3 = j + (k + 1) * a_dim1;
#line 599 "chetf2.f"
			z__4.r = d21.r * a[i__3].r - d21.i * a[i__3].i, 
				z__4.i = d21.r * a[i__3].i + d21.i * a[i__3]
				.r;
#line 599 "chetf2.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 599 "chetf2.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 599 "chetf2.f"
			wk.r = z__1.r, wk.i = z__1.i;
#line 600 "chetf2.f"
			i__2 = j + (k + 1) * a_dim1;
#line 600 "chetf2.f"
			z__3.r = d22 * a[i__2].r, z__3.i = d22 * a[i__2].i;
#line 600 "chetf2.f"
			d_cnjg(&z__5, &d21);
#line 600 "chetf2.f"
			i__3 = j + k * a_dim1;
#line 600 "chetf2.f"
			z__4.r = z__5.r * a[i__3].r - z__5.i * a[i__3].i, 
				z__4.i = z__5.r * a[i__3].i + z__5.i * a[i__3]
				.r;
#line 600 "chetf2.f"
			z__2.r = z__3.r - z__4.r, z__2.i = z__3.i - z__4.i;
#line 600 "chetf2.f"
			z__1.r = d__ * z__2.r, z__1.i = d__ * z__2.i;
#line 600 "chetf2.f"
			wkp1.r = z__1.r, wkp1.i = z__1.i;
#line 601 "chetf2.f"
			i__2 = *n;
#line 601 "chetf2.f"
			for (i__ = j; i__ <= i__2; ++i__) {
#line 602 "chetf2.f"
			    i__3 = i__ + j * a_dim1;
#line 602 "chetf2.f"
			    i__4 = i__ + j * a_dim1;
#line 602 "chetf2.f"
			    i__5 = i__ + k * a_dim1;
#line 602 "chetf2.f"
			    d_cnjg(&z__4, &wk);
#line 602 "chetf2.f"
			    z__3.r = a[i__5].r * z__4.r - a[i__5].i * z__4.i, 
				    z__3.i = a[i__5].r * z__4.i + a[i__5].i * 
				    z__4.r;
#line 602 "chetf2.f"
			    z__2.r = a[i__4].r - z__3.r, z__2.i = a[i__4].i - 
				    z__3.i;
#line 602 "chetf2.f"
			    i__6 = i__ + (k + 1) * a_dim1;
#line 602 "chetf2.f"
			    d_cnjg(&z__6, &wkp1);
#line 602 "chetf2.f"
			    z__5.r = a[i__6].r * z__6.r - a[i__6].i * z__6.i, 
				    z__5.i = a[i__6].r * z__6.i + a[i__6].i * 
				    z__6.r;
#line 602 "chetf2.f"
			    z__1.r = z__2.r - z__5.r, z__1.i = z__2.i - 
				    z__5.i;
#line 602 "chetf2.f"
			    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 604 "chetf2.f"
/* L70: */
#line 604 "chetf2.f"
			}
#line 605 "chetf2.f"
			i__2 = j + k * a_dim1;
#line 605 "chetf2.f"
			a[i__2].r = wk.r, a[i__2].i = wk.i;
#line 606 "chetf2.f"
			i__2 = j + (k + 1) * a_dim1;
#line 606 "chetf2.f"
			a[i__2].r = wkp1.r, a[i__2].i = wkp1.i;
#line 607 "chetf2.f"
			i__2 = j + j * a_dim1;
#line 607 "chetf2.f"
			i__3 = j + j * a_dim1;
#line 607 "chetf2.f"
			d__1 = a[i__3].r;
#line 607 "chetf2.f"
			z__1.r = d__1, z__1.i = 0.;
#line 607 "chetf2.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 608 "chetf2.f"
/* L80: */
#line 608 "chetf2.f"
		    }
#line 609 "chetf2.f"
		}
#line 610 "chetf2.f"
	    }
#line 611 "chetf2.f"
	}

/*        Store details of the interchanges in IPIV */

#line 615 "chetf2.f"
	if (kstep == 1) {
#line 616 "chetf2.f"
	    ipiv[k] = kp;
#line 617 "chetf2.f"
	} else {
#line 618 "chetf2.f"
	    ipiv[k] = -kp;
#line 619 "chetf2.f"
	    ipiv[k + 1] = -kp;
#line 620 "chetf2.f"
	}

/*        Increase K and return to the start of the main loop */

#line 624 "chetf2.f"
	k += kstep;
#line 625 "chetf2.f"
	goto L50;

#line 627 "chetf2.f"
    }

#line 629 "chetf2.f"
L90:
#line 630 "chetf2.f"
    return 0;

/*     End of CHETF2 */

} /* chetf2_ */

