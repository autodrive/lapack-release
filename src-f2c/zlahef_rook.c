#line 1 "zlahef_rook.f"
/* zlahef_rook.f -- translated by f2c (version 20100827).
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

#line 1 "zlahef_rook.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* \brief \b ZLAHEF_ROOK computes a partial factorization of a complex Hermitian indefinite matrix using the b
ounded Bunch-Kaufman ("rook") diagonal pivoting method (blocked algorithm, calling Level 3 BLAS). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAHEF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahef_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahef_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahef_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAHEF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KB, LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAHEF_ROOK computes a partial factorization of a complex Hermitian */
/* > matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting */
/* > method. The partial factorization has the form: */
/* > */
/* > A  =  ( I  U12 ) ( A11  0  ) (  I      0     )  if UPLO = 'U', or: */
/* >       ( 0  U22 ) (  0   D  ) ( U12**H U22**H ) */
/* > */
/* > A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L' */
/* >       ( L21  I ) (  0  A22 ) (  0      I     ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > Note that U**H denotes the conjugate transpose of U. */
/* > */
/* > ZLAHEF_ROOK is an auxiliary routine called by ZHETRF_ROOK. It uses */
/* > blocked code (calling Level 3 BLAS) to update the submatrix */
/* > A11 (if UPLO = 'U') or A22 (if UPLO = 'L'). */
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
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The maximum number of columns of the matrix A that should be */
/* >          factored.  NB should be at least 2 to allow for 2-by-2 pivot */
/* >          blocks. */
/* > \endverbatim */
/* > */
/* > \param[out] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of columns of A that were actually factored. */
/* >          KB is either NB-1 or NB, or N if N <= NB. */
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
/* >          On exit, A contains details of the partial factorization. */
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
/* >             Only the last KB elements of IPIV are set. */
/* > */
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >             interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and */
/* >             columns k and -IPIV(k) were interchanged and rows and */
/* >             columns k-1 and -IPIV(k-1) were inerchaged, */
/* >             D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
/* > */
/* >          If UPLO = 'L': */
/* >             Only the first KB elements of IPIV are set. */
/* > */
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) */
/* >             were interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and */
/* >             columns k and -IPIV(k) were interchanged and rows and */
/* >             columns k+1 and -IPIV(k+1) were inerchaged, */
/* >             D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX*16 array, dimension (LDW,NB) */
/* > \endverbatim */
/* > */
/* > \param[in] LDW */
/* > \verbatim */
/* >          LDW is INTEGER */
/* >          The leading dimension of the array W.  LDW >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization */
/* >               has been completed, but the block diagonal matrix D is */
/* >               exactly singular. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2013 */

/* > \ingroup complex16HEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2013,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zlahef_rook__(char *uplo, integer *n, integer *nb, 
	integer *kb, doublecomplex *a, integer *lda, integer *ipiv, 
	doublecomplex *w, integer *ldw, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, p;
    static doublereal t, r1;
    static doublecomplex d11, d21, d22;
    static integer jb, ii, jj, kk, kp, kw, jp1, jp2, kkw;
    static logical done;
    static integer imax, jmax;
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dtemp, sfmin;
    static integer itemp;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer kstep;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal absakk;
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal colmax;
    extern /* Subroutine */ int zlacgv_(integer *, doublecomplex *, integer *)
	    ;
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

#line 239 "zlahef_rook.f"
    /* Parameter adjustments */
#line 239 "zlahef_rook.f"
    a_dim1 = *lda;
#line 239 "zlahef_rook.f"
    a_offset = 1 + a_dim1;
#line 239 "zlahef_rook.f"
    a -= a_offset;
#line 239 "zlahef_rook.f"
    --ipiv;
#line 239 "zlahef_rook.f"
    w_dim1 = *ldw;
#line 239 "zlahef_rook.f"
    w_offset = 1 + w_dim1;
#line 239 "zlahef_rook.f"
    w -= w_offset;
#line 239 "zlahef_rook.f"

#line 239 "zlahef_rook.f"
    /* Function Body */
#line 239 "zlahef_rook.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 243 "zlahef_rook.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 247 "zlahef_rook.f"
    sfmin = dlamch_("S", (ftnlen)1);

#line 249 "zlahef_rook.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 (note that conjg(W) is actually stored) */

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

#line 257 "zlahef_rook.f"
	k = *n;
#line 258 "zlahef_rook.f"
L10:

/*        KW is the column of W which corresponds to column K of A */

#line 262 "zlahef_rook.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 266 "zlahef_rook.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 266 "zlahef_rook.f"
	    goto L30;
#line 266 "zlahef_rook.f"
	}

#line 269 "zlahef_rook.f"
	kstep = 1;
#line 270 "zlahef_rook.f"
	p = k;

/*        Copy column K of A to column KW of W and update it */

#line 274 "zlahef_rook.f"
	if (k > 1) {
#line 274 "zlahef_rook.f"
	    i__1 = k - 1;
#line 274 "zlahef_rook.f"
	    zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &
		    c__1);
#line 274 "zlahef_rook.f"
	}
#line 276 "zlahef_rook.f"
	i__1 = k + kw * w_dim1;
#line 276 "zlahef_rook.f"
	i__2 = k + k * a_dim1;
#line 276 "zlahef_rook.f"
	d__1 = a[i__2].r;
#line 276 "zlahef_rook.f"
	w[i__1].r = d__1, w[i__1].i = 0.;
#line 277 "zlahef_rook.f"
	if (k < *n) {
#line 278 "zlahef_rook.f"
	    i__1 = *n - k;
#line 278 "zlahef_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 278 "zlahef_rook.f"
	    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 280 "zlahef_rook.f"
	    i__1 = k + kw * w_dim1;
#line 280 "zlahef_rook.f"
	    i__2 = k + kw * w_dim1;
#line 280 "zlahef_rook.f"
	    d__1 = w[i__2].r;
#line 280 "zlahef_rook.f"
	    w[i__1].r = d__1, w[i__1].i = 0.;
#line 281 "zlahef_rook.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 286 "zlahef_rook.f"
	i__1 = k + kw * w_dim1;
#line 286 "zlahef_rook.f"
	absakk = (d__1 = w[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 292 "zlahef_rook.f"
	if (k > 1) {
#line 293 "zlahef_rook.f"
	    i__1 = k - 1;
#line 293 "zlahef_rook.f"
	    imax = izamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 294 "zlahef_rook.f"
	    i__1 = imax + kw * w_dim1;
#line 294 "zlahef_rook.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    kw * w_dim1]), abs(d__2));
#line 295 "zlahef_rook.f"
	} else {
#line 296 "zlahef_rook.f"
	    colmax = 0.;
#line 297 "zlahef_rook.f"
	}

#line 299 "zlahef_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 303 "zlahef_rook.f"
	    if (*info == 0) {
#line 303 "zlahef_rook.f"
		*info = k;
#line 303 "zlahef_rook.f"
	    }
#line 305 "zlahef_rook.f"
	    kp = k;
#line 306 "zlahef_rook.f"
	    i__1 = k + k * a_dim1;
#line 306 "zlahef_rook.f"
	    i__2 = k + kw * w_dim1;
#line 306 "zlahef_rook.f"
	    d__1 = w[i__2].r;
#line 306 "zlahef_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 307 "zlahef_rook.f"
	    if (k > 1) {
#line 307 "zlahef_rook.f"
		i__1 = k - 1;
#line 307 "zlahef_rook.f"
		zcopy_(&i__1, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], 
			&c__1);
#line 307 "zlahef_rook.f"
	    }
#line 309 "zlahef_rook.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */
#line 318 "zlahef_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 322 "zlahef_rook.f"
		kp = k;

#line 324 "zlahef_rook.f"
	    } else {

/*              Lop until pivot found */

#line 328 "zlahef_rook.f"
		done = FALSE_;

#line 330 "zlahef_rook.f"
L12:

/*                 BEGIN pivot search loop body */


/*                 Copy column IMAX to column KW-1 of W and update it */

#line 337 "zlahef_rook.f"
		if (imax > 1) {
#line 337 "zlahef_rook.f"
		    i__1 = imax - 1;
#line 337 "zlahef_rook.f"
		    zcopy_(&i__1, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			    w_dim1 + 1], &c__1);
#line 337 "zlahef_rook.f"
		}
#line 340 "zlahef_rook.f"
		i__1 = imax + (kw - 1) * w_dim1;
#line 340 "zlahef_rook.f"
		i__2 = imax + imax * a_dim1;
#line 340 "zlahef_rook.f"
		d__1 = a[i__2].r;
#line 340 "zlahef_rook.f"
		w[i__1].r = d__1, w[i__1].i = 0.;

#line 342 "zlahef_rook.f"
		i__1 = k - imax;
#line 342 "zlahef_rook.f"
		zcopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);
#line 344 "zlahef_rook.f"
		i__1 = k - imax;
#line 344 "zlahef_rook.f"
		zlacgv_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);

#line 346 "zlahef_rook.f"
		if (k < *n) {
#line 347 "zlahef_rook.f"
		    i__1 = *n - k;
#line 347 "zlahef_rook.f"
		    z__1.r = -1., z__1.i = -0.;
#line 347 "zlahef_rook.f"
		    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 350 "zlahef_rook.f"
		    i__1 = imax + (kw - 1) * w_dim1;
#line 350 "zlahef_rook.f"
		    i__2 = imax + (kw - 1) * w_dim1;
#line 350 "zlahef_rook.f"
		    d__1 = w[i__2].r;
#line 350 "zlahef_rook.f"
		    w[i__1].r = d__1, w[i__1].i = 0.;
#line 351 "zlahef_rook.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 357 "zlahef_rook.f"
		if (imax != k) {
#line 358 "zlahef_rook.f"
		    i__1 = k - imax;
#line 358 "zlahef_rook.f"
		    jmax = imax + izamax_(&i__1, &w[imax + 1 + (kw - 1) * 
			    w_dim1], &c__1);
#line 360 "zlahef_rook.f"
		    i__1 = jmax + (kw - 1) * w_dim1;
#line 360 "zlahef_rook.f"
		    rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    w[jmax + (kw - 1) * w_dim1]), abs(d__2));
#line 361 "zlahef_rook.f"
		} else {
#line 362 "zlahef_rook.f"
		    rowmax = 0.;
#line 363 "zlahef_rook.f"
		}

#line 365 "zlahef_rook.f"
		if (imax > 1) {
#line 366 "zlahef_rook.f"
		    i__1 = imax - 1;
#line 366 "zlahef_rook.f"
		    itemp = izamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
#line 367 "zlahef_rook.f"
		    i__1 = itemp + (kw - 1) * w_dim1;
#line 367 "zlahef_rook.f"
		    dtemp = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    itemp + (kw - 1) * w_dim1]), abs(d__2));
#line 368 "zlahef_rook.f"
		    if (dtemp > rowmax) {
#line 369 "zlahef_rook.f"
			rowmax = dtemp;
#line 370 "zlahef_rook.f"
			jmax = itemp;
#line 371 "zlahef_rook.f"
		    }
#line 372 "zlahef_rook.f"
		}

/*                 Case(2) */
/*                 Equivalent to testing for */
/*                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 379 "zlahef_rook.f"
		i__1 = imax + (kw - 1) * w_dim1;
#line 379 "zlahef_rook.f"
		if (! ((d__1 = w[i__1].r, abs(d__1)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 385 "zlahef_rook.f"
		    kp = imax;

/*                    copy column KW-1 of W to column KW of W */

#line 389 "zlahef_rook.f"
		    zcopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 391 "zlahef_rook.f"
		    done = TRUE_;

/*                 Case(3) */
/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 397 "zlahef_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 403 "zlahef_rook.f"
		    kp = imax;
#line 404 "zlahef_rook.f"
		    kstep = 2;
#line 405 "zlahef_rook.f"
		    done = TRUE_;

/*                 Case(4) */
#line 408 "zlahef_rook.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 412 "zlahef_rook.f"
		    p = imax;
#line 413 "zlahef_rook.f"
		    colmax = rowmax;
#line 414 "zlahef_rook.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 418 "zlahef_rook.f"
		    zcopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 420 "zlahef_rook.f"
		}


/*                 END pivot search loop body */

#line 425 "zlahef_rook.f"
		if (! done) {
#line 425 "zlahef_rook.f"
		    goto L12;
#line 425 "zlahef_rook.f"
		}

#line 427 "zlahef_rook.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 435 "zlahef_rook.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 439 "zlahef_rook.f"
	    kkw = *nb + kk - *n;

/*           Interchange rows and columns P and K. */
/*           Updated column P is already stored in column KW of W. */

#line 444 "zlahef_rook.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P of submatrix A */
/*              at step K. No need to copy element into columns */
/*              K and K-1 of A for 2-by-2 pivot, since these columns */
/*              will be later overwritten. */

#line 451 "zlahef_rook.f"
		i__1 = p + p * a_dim1;
#line 451 "zlahef_rook.f"
		i__2 = k + k * a_dim1;
#line 451 "zlahef_rook.f"
		d__1 = a[i__2].r;
#line 451 "zlahef_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 452 "zlahef_rook.f"
		i__1 = k - 1 - p;
#line 452 "zlahef_rook.f"
		zcopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * 
			a_dim1], lda);
#line 454 "zlahef_rook.f"
		i__1 = k - 1 - p;
#line 454 "zlahef_rook.f"
		zlacgv_(&i__1, &a[p + (p + 1) * a_dim1], lda);
#line 455 "zlahef_rook.f"
		if (p > 1) {
#line 455 "zlahef_rook.f"
		    i__1 = p - 1;
#line 455 "zlahef_rook.f"
		    zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 
			    1], &c__1);
#line 455 "zlahef_rook.f"
		}

/*              Interchange rows K and P in the last K+1 to N columns of A */
/*              (columns K and K-1 of A for 2-by-2 pivot will be */
/*              later overwritten). Interchange rows K and P */
/*              in last KKW to NB columns of W. */

#line 463 "zlahef_rook.f"
		if (k < *n) {
#line 463 "zlahef_rook.f"
		    i__1 = *n - k;
#line 463 "zlahef_rook.f"
		    zswap_(&i__1, &a[k + (k + 1) * a_dim1], lda, &a[p + (k + 
			    1) * a_dim1], lda);
#line 463 "zlahef_rook.f"
		}
#line 466 "zlahef_rook.f"
		i__1 = *n - kk + 1;
#line 466 "zlahef_rook.f"
		zswap_(&i__1, &w[k + kkw * w_dim1], ldw, &w[p + kkw * w_dim1],
			 ldw);
#line 468 "zlahef_rook.f"
	    }

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KKW of W. */

#line 473 "zlahef_rook.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K-1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 480 "zlahef_rook.f"
		i__1 = kp + kp * a_dim1;
#line 480 "zlahef_rook.f"
		i__2 = kk + kk * a_dim1;
#line 480 "zlahef_rook.f"
		d__1 = a[i__2].r;
#line 480 "zlahef_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 481 "zlahef_rook.f"
		i__1 = kk - 1 - kp;
#line 481 "zlahef_rook.f"
		zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 483 "zlahef_rook.f"
		i__1 = kk - 1 - kp;
#line 483 "zlahef_rook.f"
		zlacgv_(&i__1, &a[kp + (kp + 1) * a_dim1], lda);
#line 484 "zlahef_rook.f"
		if (kp > 1) {
#line 484 "zlahef_rook.f"
		    i__1 = kp - 1;
#line 484 "zlahef_rook.f"
		    zcopy_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 484 "zlahef_rook.f"
		}

/*              Interchange rows KK and KP in last K+1 to N columns of A */
/*              (columns K (or K and K-1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in last KKW to NB columns of W. */

#line 492 "zlahef_rook.f"
		if (k < *n) {
#line 492 "zlahef_rook.f"
		    i__1 = *n - k;
#line 492 "zlahef_rook.f"
		    zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 492 "zlahef_rook.f"
		}
#line 495 "zlahef_rook.f"
		i__1 = *n - kk + 1;
#line 495 "zlahef_rook.f"
		zswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 497 "zlahef_rook.f"
	    }

#line 499 "zlahef_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column kw of W now holds */

/*              W(kw) = U(k)*D(k), */

/*              where U(k) is the k-th column of U */

/*              (1) Store subdiag. elements of column U(k) */
/*              and 1-by-1 block D(k) in column k of A. */
/*              (NOTE: Diagonal element U(k,k) is a UNIT element */
/*              and not stored) */
/*                 A(k,k) := D(k,k) = W(k,kw) */
/*                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k) */

/*              (NOTE: No need to use for Hermitian matrix */
/*              A( K, K ) = REAL( W( K, K) ) to separately copy diagonal */
/*              element D(k,k) from W (potentially saves only one load)) */
#line 517 "zlahef_rook.f"
		zcopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 518 "zlahef_rook.f"
		if (k > 1) {

/*                 (NOTE: No need to check if A(k,k) is NOT ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  case A(k,k) = 0 falls into 2x2 pivot case(3)) */

/*                 Handle division by a small number */

#line 526 "zlahef_rook.f"
		    i__1 = k + k * a_dim1;
#line 526 "zlahef_rook.f"
		    t = a[i__1].r;
#line 527 "zlahef_rook.f"
		    if (abs(t) >= sfmin) {
#line 528 "zlahef_rook.f"
			r1 = 1. / t;
#line 529 "zlahef_rook.f"
			i__1 = k - 1;
#line 529 "zlahef_rook.f"
			zdscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 530 "zlahef_rook.f"
		    } else {
#line 531 "zlahef_rook.f"
			i__1 = k - 1;
#line 531 "zlahef_rook.f"
			for (ii = 1; ii <= i__1; ++ii) {
#line 532 "zlahef_rook.f"
			    i__2 = ii + k * a_dim1;
#line 532 "zlahef_rook.f"
			    i__3 = ii + k * a_dim1;
#line 532 "zlahef_rook.f"
			    z__1.r = a[i__3].r / t, z__1.i = a[i__3].i / t;
#line 532 "zlahef_rook.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 533 "zlahef_rook.f"
/* L14: */
#line 533 "zlahef_rook.f"
			}
#line 534 "zlahef_rook.f"
		    }

/*                 (2) Conjugate column W(kw) */

#line 538 "zlahef_rook.f"
		    i__1 = k - 1;
#line 538 "zlahef_rook.f"
		    zlacgv_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 539 "zlahef_rook.f"
		}

#line 541 "zlahef_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns kw and kw-1 of W now hold */

/*              ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              (1) Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2 */
/*              block D(k-1:k,k-1:k) in columns k-1 and k of A. */
/*              (NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT */
/*              block and not stored) */
/*                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw) */
/*                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) = */
/*                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) ) */

#line 558 "zlahef_rook.f"
		if (k > 2) {

/*                 Factor out the columns of the inverse of 2-by-2 pivot */
/*                 block D, so that each column contains 1, to reduce the */
/*                 number of FLOPS when we multiply panel */
/*                 ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1). */

/*                 D**(-1) = ( d11 cj(d21) )**(-1) = */
/*                           ( d21    d22 ) */

/*                 = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) = */
/*                                          ( (-d21) (     d11 ) ) */

/*                 = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) * */

/*                   * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) = */
/*                     (     (      -1 )           ( d11/conj(d21) ) ) */

/*                 = 1/(|d21|**2) * 1/(D22*D11-1) * */

/*                   * ( d21*( D11 ) conj(d21)*(  -1 ) ) = */
/*                     (     (  -1 )           ( D22 ) ) */

/*                 = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) = */
/*                                      (     (  -1 )           ( D22 ) ) */

/*                 = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) = */
/*                   (               (  -1 )         ( D22 ) ) */

/*                 Handle division by a small number. (NOTE: order of */
/*                 operations is important) */

/*                 = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) ) */
/*                   (   ((  -1 )          )   (( D22 )     ) ), */

/*                 where D11 = d22/d21, */
/*                       D22 = d11/conj(d21), */
/*                       D21 = d21, */
/*                       T = 1/(D22*D11-1). */

/*                 (NOTE: No need to check for division by ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  (a) d21 != 0 in 2x2 pivot case(4), */
/*                      since |d21| should be larger than |d11| and |d22|; */
/*                  (b) (D22*D11 - 1) != 0, since from (a), */
/*                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.) */

#line 605 "zlahef_rook.f"
		    i__1 = k - 1 + kw * w_dim1;
#line 605 "zlahef_rook.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 606 "zlahef_rook.f"
		    d_cnjg(&z__2, &d21);
#line 606 "zlahef_rook.f"
		    z_div(&z__1, &w[k + kw * w_dim1], &z__2);
#line 606 "zlahef_rook.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 607 "zlahef_rook.f"
		    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d21);
#line 607 "zlahef_rook.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 608 "zlahef_rook.f"
		    z__1.r = d11.r * d22.r - d11.i * d22.i, z__1.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 608 "zlahef_rook.f"
		    t = 1. / (z__1.r - 1.);

/*                 Update elements in columns A(k-1) and A(k) as */
/*                 dot products of rows of ( W(kw-1) W(kw) ) and columns */
/*                 of D**(-1) */

#line 614 "zlahef_rook.f"
		    i__1 = k - 2;
#line 614 "zlahef_rook.f"
		    for (j = 1; j <= i__1; ++j) {
#line 615 "zlahef_rook.f"
			i__2 = j + (k - 1) * a_dim1;
#line 615 "zlahef_rook.f"
			i__3 = j + (kw - 1) * w_dim1;
#line 615 "zlahef_rook.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 615 "zlahef_rook.f"
			i__4 = j + kw * w_dim1;
#line 615 "zlahef_rook.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 615 "zlahef_rook.f"
			z_div(&z__2, &z__3, &d21);
#line 615 "zlahef_rook.f"
			z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 615 "zlahef_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 617 "zlahef_rook.f"
			i__2 = j + k * a_dim1;
#line 617 "zlahef_rook.f"
			i__3 = j + kw * w_dim1;
#line 617 "zlahef_rook.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 617 "zlahef_rook.f"
			i__4 = j + (kw - 1) * w_dim1;
#line 617 "zlahef_rook.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 617 "zlahef_rook.f"
			d_cnjg(&z__5, &d21);
#line 617 "zlahef_rook.f"
			z_div(&z__2, &z__3, &z__5);
#line 617 "zlahef_rook.f"
			z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 617 "zlahef_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 619 "zlahef_rook.f"
/* L20: */
#line 619 "zlahef_rook.f"
		    }
#line 620 "zlahef_rook.f"
		}

/*              Copy D(k) to A */

#line 624 "zlahef_rook.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 624 "zlahef_rook.f"
		i__2 = k - 1 + (kw - 1) * w_dim1;
#line 624 "zlahef_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 625 "zlahef_rook.f"
		i__1 = k - 1 + k * a_dim1;
#line 625 "zlahef_rook.f"
		i__2 = k - 1 + kw * w_dim1;
#line 625 "zlahef_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 626 "zlahef_rook.f"
		i__1 = k + k * a_dim1;
#line 626 "zlahef_rook.f"
		i__2 = k + kw * w_dim1;
#line 626 "zlahef_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

/*              (2) Conjugate columns W(kw) and W(kw-1) */

#line 630 "zlahef_rook.f"
		i__1 = k - 1;
#line 630 "zlahef_rook.f"
		zlacgv_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 631 "zlahef_rook.f"
		i__1 = k - 2;
#line 631 "zlahef_rook.f"
		zlacgv_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);

#line 633 "zlahef_rook.f"
	    }

#line 635 "zlahef_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 639 "zlahef_rook.f"
	if (kstep == 1) {
#line 640 "zlahef_rook.f"
	    ipiv[k] = kp;
#line 641 "zlahef_rook.f"
	} else {
#line 642 "zlahef_rook.f"
	    ipiv[k] = -p;
#line 643 "zlahef_rook.f"
	    ipiv[k - 1] = -kp;
#line 644 "zlahef_rook.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 648 "zlahef_rook.f"
	k -= kstep;
#line 649 "zlahef_rook.f"
	goto L10;

#line 651 "zlahef_rook.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**H = A11 - U12*W**H */

/*        computing blocks of NB columns at a time (note that conjg(W) is */
/*        actually stored) */

#line 660 "zlahef_rook.f"
	i__1 = -(*nb);
#line 660 "zlahef_rook.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 661 "zlahef_rook.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 661 "zlahef_rook.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 665 "zlahef_rook.f"
	    i__2 = j + jb - 1;
#line 665 "zlahef_rook.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 666 "zlahef_rook.f"
		i__3 = jj + jj * a_dim1;
#line 666 "zlahef_rook.f"
		i__4 = jj + jj * a_dim1;
#line 666 "zlahef_rook.f"
		d__1 = a[i__4].r;
#line 666 "zlahef_rook.f"
		a[i__3].r = d__1, a[i__3].i = 0.;
#line 667 "zlahef_rook.f"
		i__3 = jj - j + 1;
#line 667 "zlahef_rook.f"
		i__4 = *n - k;
#line 667 "zlahef_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 667 "zlahef_rook.f"
		zgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, 
			&a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 670 "zlahef_rook.f"
		i__3 = jj + jj * a_dim1;
#line 670 "zlahef_rook.f"
		i__4 = jj + jj * a_dim1;
#line 670 "zlahef_rook.f"
		d__1 = a[i__4].r;
#line 670 "zlahef_rook.f"
		a[i__3].r = d__1, a[i__3].i = 0.;
#line 671 "zlahef_rook.f"
/* L40: */
#line 671 "zlahef_rook.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 675 "zlahef_rook.f"
	    if (j >= 2) {
#line 675 "zlahef_rook.f"
		i__2 = j - 1;
#line 675 "zlahef_rook.f"
		i__3 = *n - k;
#line 675 "zlahef_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 675 "zlahef_rook.f"
		zgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, 
			&a[(k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * 
			w_dim1], ldw, &c_b1, &a[j * a_dim1 + 1], lda, (ftnlen)
			12, (ftnlen)9);
#line 675 "zlahef_rook.f"
	    }
#line 679 "zlahef_rook.f"
/* L50: */
#line 679 "zlahef_rook.f"
	}

/*        Put U12 in standard form by partially undoing the interchanges */
/*        in of rows in columns k+1:n looping backwards from k+1 to n */

#line 684 "zlahef_rook.f"
	j = k + 1;
#line 685 "zlahef_rook.f"
L60:

/*           Undo the interchanges (if any) of rows J and JP2 */
/*           (or J and JP2, and J+1 and JP1) at each step J */

#line 690 "zlahef_rook.f"
	kstep = 1;
#line 691 "zlahef_rook.f"
	jp1 = 1;
/*           (Here, J is a diagonal index) */
#line 693 "zlahef_rook.f"
	jj = j;
#line 694 "zlahef_rook.f"
	jp2 = ipiv[j];
#line 695 "zlahef_rook.f"
	if (jp2 < 0) {
#line 696 "zlahef_rook.f"
	    jp2 = -jp2;
/*              (Here, J is a diagonal index) */
#line 698 "zlahef_rook.f"
	    ++j;
#line 699 "zlahef_rook.f"
	    jp1 = -ipiv[j];
#line 700 "zlahef_rook.f"
	    kstep = 2;
#line 701 "zlahef_rook.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length N-J+1 */
/*           of the rows to swap back doesn't include diagonal element) */
#line 704 "zlahef_rook.f"
	++j;
#line 705 "zlahef_rook.f"
	if (jp2 != jj && j <= *n) {
#line 705 "zlahef_rook.f"
	    i__1 = *n - j + 1;
#line 705 "zlahef_rook.f"
	    zswap_(&i__1, &a[jp2 + j * a_dim1], lda, &a[jj + j * a_dim1], lda)
		    ;
#line 705 "zlahef_rook.f"
	}
#line 707 "zlahef_rook.f"
	++jj;
#line 708 "zlahef_rook.f"
	if (kstep == 2 && jp1 != jj && j <= *n) {
#line 708 "zlahef_rook.f"
	    i__1 = *n - j + 1;
#line 708 "zlahef_rook.f"
	    zswap_(&i__1, &a[jp1 + j * a_dim1], lda, &a[jj + j * a_dim1], lda)
		    ;
#line 708 "zlahef_rook.f"
	}
#line 710 "zlahef_rook.f"
	if (j < *n) {
#line 710 "zlahef_rook.f"
	    goto L60;
#line 710 "zlahef_rook.f"
	}

/*        Set KB to the number of columns factorized */

#line 715 "zlahef_rook.f"
	*kb = *n - k;

#line 717 "zlahef_rook.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 (note that conjg(W) is actually stored) */

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 725 "zlahef_rook.f"
	k = 1;
#line 726 "zlahef_rook.f"
L70:

/*        Exit from loop */

#line 730 "zlahef_rook.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 730 "zlahef_rook.f"
	    goto L90;
#line 730 "zlahef_rook.f"
	}

#line 733 "zlahef_rook.f"
	kstep = 1;
#line 734 "zlahef_rook.f"
	p = k;

/*        Copy column K of A to column K of W and update column K of W */

#line 738 "zlahef_rook.f"
	i__1 = k + k * w_dim1;
#line 738 "zlahef_rook.f"
	i__2 = k + k * a_dim1;
#line 738 "zlahef_rook.f"
	d__1 = a[i__2].r;
#line 738 "zlahef_rook.f"
	w[i__1].r = d__1, w[i__1].i = 0.;
#line 739 "zlahef_rook.f"
	if (k < *n) {
#line 739 "zlahef_rook.f"
	    i__1 = *n - k;
#line 739 "zlahef_rook.f"
	    zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &w[k + 1 + k * 
		    w_dim1], &c__1);
#line 739 "zlahef_rook.f"
	}
#line 741 "zlahef_rook.f"
	if (k > 1) {
#line 742 "zlahef_rook.f"
	    i__1 = *n - k + 1;
#line 742 "zlahef_rook.f"
	    i__2 = k - 1;
#line 742 "zlahef_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 742 "zlahef_rook.f"
	    zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &
		    w[k + w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1, (
		    ftnlen)12);
#line 744 "zlahef_rook.f"
	    i__1 = k + k * w_dim1;
#line 744 "zlahef_rook.f"
	    i__2 = k + k * w_dim1;
#line 744 "zlahef_rook.f"
	    d__1 = w[i__2].r;
#line 744 "zlahef_rook.f"
	    w[i__1].r = d__1, w[i__1].i = 0.;
#line 745 "zlahef_rook.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 750 "zlahef_rook.f"
	i__1 = k + k * w_dim1;
#line 750 "zlahef_rook.f"
	absakk = (d__1 = w[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 756 "zlahef_rook.f"
	if (k < *n) {
#line 757 "zlahef_rook.f"
	    i__1 = *n - k;
#line 757 "zlahef_rook.f"
	    imax = k + izamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 758 "zlahef_rook.f"
	    i__1 = imax + k * w_dim1;
#line 758 "zlahef_rook.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    k * w_dim1]), abs(d__2));
#line 759 "zlahef_rook.f"
	} else {
#line 760 "zlahef_rook.f"
	    colmax = 0.;
#line 761 "zlahef_rook.f"
	}

#line 763 "zlahef_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 767 "zlahef_rook.f"
	    if (*info == 0) {
#line 767 "zlahef_rook.f"
		*info = k;
#line 767 "zlahef_rook.f"
	    }
#line 769 "zlahef_rook.f"
	    kp = k;
#line 770 "zlahef_rook.f"
	    i__1 = k + k * a_dim1;
#line 770 "zlahef_rook.f"
	    i__2 = k + k * w_dim1;
#line 770 "zlahef_rook.f"
	    d__1 = w[i__2].r;
#line 770 "zlahef_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 771 "zlahef_rook.f"
	    if (k < *n) {
#line 771 "zlahef_rook.f"
		i__1 = *n - k;
#line 771 "zlahef_rook.f"
		zcopy_(&i__1, &w[k + 1 + k * w_dim1], &c__1, &a[k + 1 + k * 
			a_dim1], &c__1);
#line 771 "zlahef_rook.f"
	    }
#line 773 "zlahef_rook.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 783 "zlahef_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 787 "zlahef_rook.f"
		kp = k;

#line 789 "zlahef_rook.f"
	    } else {

#line 791 "zlahef_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 795 "zlahef_rook.f"
L72:

/*                 BEGIN pivot search loop body */


/*                 Copy column IMAX to column k+1 of W and update it */

#line 802 "zlahef_rook.f"
		i__1 = imax - k;
#line 802 "zlahef_rook.f"
		zcopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 803 "zlahef_rook.f"
		i__1 = imax - k;
#line 803 "zlahef_rook.f"
		zlacgv_(&i__1, &w[k + (k + 1) * w_dim1], &c__1);
#line 804 "zlahef_rook.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 804 "zlahef_rook.f"
		i__2 = imax + imax * a_dim1;
#line 804 "zlahef_rook.f"
		d__1 = a[i__2].r;
#line 804 "zlahef_rook.f"
		w[i__1].r = d__1, w[i__1].i = 0.;

#line 806 "zlahef_rook.f"
		if (imax < *n) {
#line 806 "zlahef_rook.f"
		    i__1 = *n - imax;
#line 806 "zlahef_rook.f"
		    zcopy_(&i__1, &a[imax + 1 + imax * a_dim1], &c__1, &w[
			    imax + 1 + (k + 1) * w_dim1], &c__1);
#line 806 "zlahef_rook.f"
		}

#line 810 "zlahef_rook.f"
		if (k > 1) {
#line 811 "zlahef_rook.f"
		    i__1 = *n - k + 1;
#line 811 "zlahef_rook.f"
		    i__2 = k - 1;
#line 811 "zlahef_rook.f"
		    z__1.r = -1., z__1.i = -0.;
#line 811 "zlahef_rook.f"
		    zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1]
			    , lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 
			    1) * w_dim1], &c__1, (ftnlen)12);
#line 814 "zlahef_rook.f"
		    i__1 = imax + (k + 1) * w_dim1;
#line 814 "zlahef_rook.f"
		    i__2 = imax + (k + 1) * w_dim1;
#line 814 "zlahef_rook.f"
		    d__1 = w[i__2].r;
#line 814 "zlahef_rook.f"
		    w[i__1].r = d__1, w[i__1].i = 0.;
#line 815 "zlahef_rook.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 821 "zlahef_rook.f"
		if (imax != k) {
#line 822 "zlahef_rook.f"
		    i__1 = imax - k;
#line 822 "zlahef_rook.f"
		    jmax = k - 1 + izamax_(&i__1, &w[k + (k + 1) * w_dim1], &
			    c__1);
#line 823 "zlahef_rook.f"
		    i__1 = jmax + (k + 1) * w_dim1;
#line 823 "zlahef_rook.f"
		    rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    w[jmax + (k + 1) * w_dim1]), abs(d__2));
#line 824 "zlahef_rook.f"
		} else {
#line 825 "zlahef_rook.f"
		    rowmax = 0.;
#line 826 "zlahef_rook.f"
		}

#line 828 "zlahef_rook.f"
		if (imax < *n) {
#line 829 "zlahef_rook.f"
		    i__1 = *n - imax;
#line 829 "zlahef_rook.f"
		    itemp = imax + izamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
#line 830 "zlahef_rook.f"
		    i__1 = itemp + (k + 1) * w_dim1;
#line 830 "zlahef_rook.f"
		    dtemp = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    itemp + (k + 1) * w_dim1]), abs(d__2));
#line 831 "zlahef_rook.f"
		    if (dtemp > rowmax) {
#line 832 "zlahef_rook.f"
			rowmax = dtemp;
#line 833 "zlahef_rook.f"
			jmax = itemp;
#line 834 "zlahef_rook.f"
		    }
#line 835 "zlahef_rook.f"
		}

/*                 Case(2) */
/*                 Equivalent to testing for */
/*                 ABS( REAL( W( IMAX,K+1 ) ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 842 "zlahef_rook.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 842 "zlahef_rook.f"
		if (! ((d__1 = w[i__1].r, abs(d__1)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 848 "zlahef_rook.f"
		    kp = imax;

/*                    copy column K+1 of W to column K of W */

#line 852 "zlahef_rook.f"
		    i__1 = *n - k + 1;
#line 852 "zlahef_rook.f"
		    zcopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 854 "zlahef_rook.f"
		    done = TRUE_;

/*                 Case(3) */
/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 860 "zlahef_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 866 "zlahef_rook.f"
		    kp = imax;
#line 867 "zlahef_rook.f"
		    kstep = 2;
#line 868 "zlahef_rook.f"
		    done = TRUE_;

/*                 Case(4) */
#line 871 "zlahef_rook.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 875 "zlahef_rook.f"
		    p = imax;
#line 876 "zlahef_rook.f"
		    colmax = rowmax;
#line 877 "zlahef_rook.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 881 "zlahef_rook.f"
		    i__1 = *n - k + 1;
#line 881 "zlahef_rook.f"
		    zcopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 883 "zlahef_rook.f"
		}


/*                 End pivot search loop body */

#line 888 "zlahef_rook.f"
		if (! done) {
#line 888 "zlahef_rook.f"
		    goto L72;
#line 888 "zlahef_rook.f"
		}

#line 890 "zlahef_rook.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 898 "zlahef_rook.f"
	    kk = k + kstep - 1;

/*           Interchange rows and columns P and K (only for 2-by-2 pivot). */
/*           Updated column P is already stored in column K of W. */

#line 903 "zlahef_rook.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column KK-1 to column P of submatrix A */
/*              at step K. No need to copy element into columns */
/*              K and K+1 of A for 2-by-2 pivot, since these columns */
/*              will be later overwritten. */

#line 910 "zlahef_rook.f"
		i__1 = p + p * a_dim1;
#line 910 "zlahef_rook.f"
		i__2 = k + k * a_dim1;
#line 910 "zlahef_rook.f"
		d__1 = a[i__2].r;
#line 910 "zlahef_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 911 "zlahef_rook.f"
		i__1 = p - k - 1;
#line 911 "zlahef_rook.f"
		zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[p + (k + 1) * 
			a_dim1], lda);
#line 912 "zlahef_rook.f"
		i__1 = p - k - 1;
#line 912 "zlahef_rook.f"
		zlacgv_(&i__1, &a[p + (k + 1) * a_dim1], lda);
#line 913 "zlahef_rook.f"
		if (p < *n) {
#line 913 "zlahef_rook.f"
		    i__1 = *n - p;
#line 913 "zlahef_rook.f"
		    zcopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p 
			    * a_dim1], &c__1);
#line 913 "zlahef_rook.f"
		}

/*              Interchange rows K and P in first K-1 columns of A */
/*              (columns K and K+1 of A for 2-by-2 pivot will be */
/*              later overwritten). Interchange rows K and P */
/*              in first KK columns of W. */

#line 921 "zlahef_rook.f"
		if (k > 1) {
#line 921 "zlahef_rook.f"
		    i__1 = k - 1;
#line 921 "zlahef_rook.f"
		    zswap_(&i__1, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 921 "zlahef_rook.f"
		}
#line 923 "zlahef_rook.f"
		zswap_(&kk, &w[k + w_dim1], ldw, &w[p + w_dim1], ldw);
#line 924 "zlahef_rook.f"
	    }

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KK of W. */

#line 929 "zlahef_rook.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K+1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 936 "zlahef_rook.f"
		i__1 = kp + kp * a_dim1;
#line 936 "zlahef_rook.f"
		i__2 = kk + kk * a_dim1;
#line 936 "zlahef_rook.f"
		d__1 = a[i__2].r;
#line 936 "zlahef_rook.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 937 "zlahef_rook.f"
		i__1 = kp - kk - 1;
#line 937 "zlahef_rook.f"
		zcopy_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 939 "zlahef_rook.f"
		i__1 = kp - kk - 1;
#line 939 "zlahef_rook.f"
		zlacgv_(&i__1, &a[kp + (kk + 1) * a_dim1], lda);
#line 940 "zlahef_rook.f"
		if (kp < *n) {
#line 940 "zlahef_rook.f"
		    i__1 = *n - kp;
#line 940 "zlahef_rook.f"
		    zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 940 "zlahef_rook.f"
		}

/*              Interchange rows KK and KP in first K-1 columns of A */
/*              (column K (or K and K+1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in first KK columns of W. */

#line 948 "zlahef_rook.f"
		if (k > 1) {
#line 948 "zlahef_rook.f"
		    i__1 = k - 1;
#line 948 "zlahef_rook.f"
		    zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 948 "zlahef_rook.f"
		}
#line 950 "zlahef_rook.f"
		zswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 951 "zlahef_rook.f"
	    }

#line 953 "zlahef_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of W now holds */

/*              W(k) = L(k)*D(k), */

/*              where L(k) is the k-th column of L */

/*              (1) Store subdiag. elements of column L(k) */
/*              and 1-by-1 block D(k) in column k of A. */
/*              (NOTE: Diagonal element L(k,k) is a UNIT element */
/*              and not stored) */
/*                 A(k,k) := D(k,k) = W(k,k) */
/*                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k) */

/*              (NOTE: No need to use for Hermitian matrix */
/*              A( K, K ) = REAL( W( K, K) ) to separately copy diagonal */
/*              element D(k,k) from W (potentially saves only one load)) */
#line 971 "zlahef_rook.f"
		i__1 = *n - k + 1;
#line 971 "zlahef_rook.f"
		zcopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 972 "zlahef_rook.f"
		if (k < *n) {

/*                 (NOTE: No need to check if A(k,k) is NOT ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  case A(k,k) = 0 falls into 2x2 pivot case(3)) */

/*                 Handle division by a small number */

#line 980 "zlahef_rook.f"
		    i__1 = k + k * a_dim1;
#line 980 "zlahef_rook.f"
		    t = a[i__1].r;
#line 981 "zlahef_rook.f"
		    if (abs(t) >= sfmin) {
#line 982 "zlahef_rook.f"
			r1 = 1. / t;
#line 983 "zlahef_rook.f"
			i__1 = *n - k;
#line 983 "zlahef_rook.f"
			zdscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 984 "zlahef_rook.f"
		    } else {
#line 985 "zlahef_rook.f"
			i__1 = *n;
#line 985 "zlahef_rook.f"
			for (ii = k + 1; ii <= i__1; ++ii) {
#line 986 "zlahef_rook.f"
			    i__2 = ii + k * a_dim1;
#line 986 "zlahef_rook.f"
			    i__3 = ii + k * a_dim1;
#line 986 "zlahef_rook.f"
			    z__1.r = a[i__3].r / t, z__1.i = a[i__3].i / t;
#line 986 "zlahef_rook.f"
			    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 987 "zlahef_rook.f"
/* L74: */
#line 987 "zlahef_rook.f"
			}
#line 988 "zlahef_rook.f"
		    }

/*                 (2) Conjugate column W(k) */

#line 992 "zlahef_rook.f"
		    i__1 = *n - k;
#line 992 "zlahef_rook.f"
		    zlacgv_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 993 "zlahef_rook.f"
		}

#line 995 "zlahef_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

/*              (1) Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2 */
/*              block D(k:k+1,k:k+1) in columns k and k+1 of A. */
/*              NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT */
/*              block and not stored. */
/*                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1) */
/*                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) = */
/*                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) ) */

#line 1012 "zlahef_rook.f"
		if (k < *n - 1) {

/*                 Factor out the columns of the inverse of 2-by-2 pivot */
/*                 block D, so that each column contains 1, to reduce the */
/*                 number of FLOPS when we multiply panel */
/*                 ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1). */

/*                 D**(-1) = ( d11 cj(d21) )**(-1) = */
/*                           ( d21    d22 ) */

/*                 = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) = */
/*                                          ( (-d21) (     d11 ) ) */

/*                 = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) * */

/*                   * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) = */
/*                     (     (      -1 )           ( d11/conj(d21) ) ) */

/*                 = 1/(|d21|**2) * 1/(D22*D11-1) * */

/*                   * ( d21*( D11 ) conj(d21)*(  -1 ) ) = */
/*                     (     (  -1 )           ( D22 ) ) */

/*                 = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) = */
/*                                      (     (  -1 )           ( D22 ) ) */

/*                 = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) = */
/*                   (               (  -1 )         ( D22 ) ) */

/*                 Handle division by a small number. (NOTE: order of */
/*                 operations is important) */

/*                 = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) ) */
/*                   (   ((  -1 )          )   (( D22 )     ) ), */

/*                 where D11 = d22/d21, */
/*                       D22 = d11/conj(d21), */
/*                       D21 = d21, */
/*                       T = 1/(D22*D11-1). */

/*                 (NOTE: No need to check for division by ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  (a) d21 != 0 in 2x2 pivot case(4), */
/*                      since |d21| should be larger than |d11| and |d22|; */
/*                  (b) (D22*D11 - 1) != 0, since from (a), */
/*                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.) */

#line 1059 "zlahef_rook.f"
		    i__1 = k + 1 + k * w_dim1;
#line 1059 "zlahef_rook.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 1060 "zlahef_rook.f"
		    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
#line 1060 "zlahef_rook.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 1061 "zlahef_rook.f"
		    d_cnjg(&z__2, &d21);
#line 1061 "zlahef_rook.f"
		    z_div(&z__1, &w[k + k * w_dim1], &z__2);
#line 1061 "zlahef_rook.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 1062 "zlahef_rook.f"
		    z__1.r = d11.r * d22.r - d11.i * d22.i, z__1.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 1062 "zlahef_rook.f"
		    t = 1. / (z__1.r - 1.);

/*                 Update elements in columns A(k) and A(k+1) as */
/*                 dot products of rows of ( W(k) W(k+1) ) and columns */
/*                 of D**(-1) */

#line 1068 "zlahef_rook.f"
		    i__1 = *n;
#line 1068 "zlahef_rook.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 1069 "zlahef_rook.f"
			i__2 = j + k * a_dim1;
#line 1069 "zlahef_rook.f"
			i__3 = j + k * w_dim1;
#line 1069 "zlahef_rook.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 1069 "zlahef_rook.f"
			i__4 = j + (k + 1) * w_dim1;
#line 1069 "zlahef_rook.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 1069 "zlahef_rook.f"
			d_cnjg(&z__5, &d21);
#line 1069 "zlahef_rook.f"
			z_div(&z__2, &z__3, &z__5);
#line 1069 "zlahef_rook.f"
			z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 1069 "zlahef_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 1071 "zlahef_rook.f"
			i__2 = j + (k + 1) * a_dim1;
#line 1071 "zlahef_rook.f"
			i__3 = j + (k + 1) * w_dim1;
#line 1071 "zlahef_rook.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 1071 "zlahef_rook.f"
			i__4 = j + k * w_dim1;
#line 1071 "zlahef_rook.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 1071 "zlahef_rook.f"
			z_div(&z__2, &z__3, &d21);
#line 1071 "zlahef_rook.f"
			z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 1071 "zlahef_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 1073 "zlahef_rook.f"
/* L80: */
#line 1073 "zlahef_rook.f"
		    }
#line 1074 "zlahef_rook.f"
		}

/*              Copy D(k) to A */

#line 1078 "zlahef_rook.f"
		i__1 = k + k * a_dim1;
#line 1078 "zlahef_rook.f"
		i__2 = k + k * w_dim1;
#line 1078 "zlahef_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 1079 "zlahef_rook.f"
		i__1 = k + 1 + k * a_dim1;
#line 1079 "zlahef_rook.f"
		i__2 = k + 1 + k * w_dim1;
#line 1079 "zlahef_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 1080 "zlahef_rook.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 1080 "zlahef_rook.f"
		i__2 = k + 1 + (k + 1) * w_dim1;
#line 1080 "zlahef_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

/*              (2) Conjugate columns W(k) and W(k+1) */

#line 1084 "zlahef_rook.f"
		i__1 = *n - k;
#line 1084 "zlahef_rook.f"
		zlacgv_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 1085 "zlahef_rook.f"
		i__1 = *n - k - 1;
#line 1085 "zlahef_rook.f"
		zlacgv_(&i__1, &w[k + 2 + (k + 1) * w_dim1], &c__1);

#line 1087 "zlahef_rook.f"
	    }

#line 1089 "zlahef_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 1093 "zlahef_rook.f"
	if (kstep == 1) {
#line 1094 "zlahef_rook.f"
	    ipiv[k] = kp;
#line 1095 "zlahef_rook.f"
	} else {
#line 1096 "zlahef_rook.f"
	    ipiv[k] = -p;
#line 1097 "zlahef_rook.f"
	    ipiv[k + 1] = -kp;
#line 1098 "zlahef_rook.f"
	}

/*        Increase K and return to the start of the main loop */

#line 1102 "zlahef_rook.f"
	k += kstep;
#line 1103 "zlahef_rook.f"
	goto L70;

#line 1105 "zlahef_rook.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**H = A22 - L21*W**H */

/*        computing blocks of NB columns at a time (note that conjg(W) is */
/*        actually stored) */

#line 1114 "zlahef_rook.f"
	i__1 = *n;
#line 1114 "zlahef_rook.f"
	i__2 = *nb;
#line 1114 "zlahef_rook.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 1115 "zlahef_rook.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 1115 "zlahef_rook.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 1119 "zlahef_rook.f"
	    i__3 = j + jb - 1;
#line 1119 "zlahef_rook.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 1120 "zlahef_rook.f"
		i__4 = jj + jj * a_dim1;
#line 1120 "zlahef_rook.f"
		i__5 = jj + jj * a_dim1;
#line 1120 "zlahef_rook.f"
		d__1 = a[i__5].r;
#line 1120 "zlahef_rook.f"
		a[i__4].r = d__1, a[i__4].i = 0.;
#line 1121 "zlahef_rook.f"
		i__4 = j + jb - jj;
#line 1121 "zlahef_rook.f"
		i__5 = k - 1;
#line 1121 "zlahef_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 1121 "zlahef_rook.f"
		zgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1]
			, &c__1, (ftnlen)12);
#line 1124 "zlahef_rook.f"
		i__4 = jj + jj * a_dim1;
#line 1124 "zlahef_rook.f"
		i__5 = jj + jj * a_dim1;
#line 1124 "zlahef_rook.f"
		d__1 = a[i__5].r;
#line 1124 "zlahef_rook.f"
		a[i__4].r = d__1, a[i__4].i = 0.;
#line 1125 "zlahef_rook.f"
/* L100: */
#line 1125 "zlahef_rook.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 1129 "zlahef_rook.f"
	    if (j + jb <= *n) {
#line 1129 "zlahef_rook.f"
		i__3 = *n - j - jb + 1;
#line 1129 "zlahef_rook.f"
		i__4 = k - 1;
#line 1129 "zlahef_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 1129 "zlahef_rook.f"
		zgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, 
			&a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 1129 "zlahef_rook.f"
	    }
#line 1133 "zlahef_rook.f"
/* L110: */
#line 1133 "zlahef_rook.f"
	}

/*        Put L21 in standard form by partially undoing the interchanges */
/*        of rows in columns 1:k-1 looping backwards from k-1 to 1 */

#line 1138 "zlahef_rook.f"
	j = k - 1;
#line 1139 "zlahef_rook.f"
L120:

/*           Undo the interchanges (if any) of rows J and JP2 */
/*           (or J and JP2, and J-1 and JP1) at each step J */

#line 1144 "zlahef_rook.f"
	kstep = 1;
#line 1145 "zlahef_rook.f"
	jp1 = 1;
/*           (Here, J is a diagonal index) */
#line 1147 "zlahef_rook.f"
	jj = j;
#line 1148 "zlahef_rook.f"
	jp2 = ipiv[j];
#line 1149 "zlahef_rook.f"
	if (jp2 < 0) {
#line 1150 "zlahef_rook.f"
	    jp2 = -jp2;
/*              (Here, J is a diagonal index) */
#line 1152 "zlahef_rook.f"
	    --j;
#line 1153 "zlahef_rook.f"
	    jp1 = -ipiv[j];
#line 1154 "zlahef_rook.f"
	    kstep = 2;
#line 1155 "zlahef_rook.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length J */
/*           of the rows to swap back doesn't include diagonal element) */
#line 1158 "zlahef_rook.f"
	--j;
#line 1159 "zlahef_rook.f"
	if (jp2 != jj && j >= 1) {
#line 1159 "zlahef_rook.f"
	    zswap_(&j, &a[jp2 + a_dim1], lda, &a[jj + a_dim1], lda);
#line 1159 "zlahef_rook.f"
	}
#line 1161 "zlahef_rook.f"
	--jj;
#line 1162 "zlahef_rook.f"
	if (kstep == 2 && jp1 != jj && j >= 1) {
#line 1162 "zlahef_rook.f"
	    zswap_(&j, &a[jp1 + a_dim1], lda, &a[jj + a_dim1], lda);
#line 1162 "zlahef_rook.f"
	}
#line 1164 "zlahef_rook.f"
	if (j > 1) {
#line 1164 "zlahef_rook.f"
	    goto L120;
#line 1164 "zlahef_rook.f"
	}

/*        Set KB to the number of columns factorized */

#line 1169 "zlahef_rook.f"
	*kb = k - 1;

#line 1171 "zlahef_rook.f"
    }
#line 1172 "zlahef_rook.f"
    return 0;

/*     End of ZLAHEF_ROOK */

} /* zlahef_rook__ */

