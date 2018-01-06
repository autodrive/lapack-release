#line 1 "clahef.f"
/* clahef.f -- translated by f2c (version 20100827).
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

#line 1 "clahef.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLAHEF computes a partial factorization of a complex Hermitian indefinite matrix using the Bunc
h-Kaufman diagonal pivoting method (blocked algorithm, calling Level 3 BLAS). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAHEF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clahef.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clahef.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clahef.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAHEF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KB, LDA, LDW, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), W( LDW, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAHEF computes a partial factorization of a complex Hermitian */
/* > matrix A using the Bunch-Kaufman diagonal pivoting method. The */
/* > partial factorization has the form: */
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
/* > CLAHEF is an auxiliary routine called by CHETRF. It uses blocked code */
/* > (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or */
/* > A22 (if UPLO = 'L'). */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >             If IPIV(k) = IPIV(k-1) < 0, then rows and columns */
/* >             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* >             is a 2-by-2 diagonal block. */
/* > */
/* >          If UPLO = 'L': */
/* >             Only the first KB elements of IPIV are set. */
/* > */
/* >             If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >             interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* >             If IPIV(k) = IPIV(k+1) < 0, then rows and columns */
/* >             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1) */
/* >             is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is COMPLEX array, dimension (LDW,NB) */
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

/* > \ingroup complexHEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2013,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int clahef_(char *uplo, integer *n, integer *nb, integer *kb,
	 doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, 
	integer *ldw, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublereal t, r1;
    static doublecomplex d11, d21, d22;
    static integer jb, jj, kk, jp, kp, kw, kkw, imax, jmax;
    static doublereal alpha;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ccopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer kstep;
    static doublereal absakk;
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    ;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal colmax, rowmax;


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

#line 229 "clahef.f"
    /* Parameter adjustments */
#line 229 "clahef.f"
    a_dim1 = *lda;
#line 229 "clahef.f"
    a_offset = 1 + a_dim1;
#line 229 "clahef.f"
    a -= a_offset;
#line 229 "clahef.f"
    --ipiv;
#line 229 "clahef.f"
    w_dim1 = *ldw;
#line 229 "clahef.f"
    w_offset = 1 + w_dim1;
#line 229 "clahef.f"
    w -= w_offset;
#line 229 "clahef.f"

#line 229 "clahef.f"
    /* Function Body */
#line 229 "clahef.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 233 "clahef.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 235 "clahef.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 (note that conjg(W) is actually stored) */

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

#line 243 "clahef.f"
	k = *n;
#line 244 "clahef.f"
L10:

/*        KW is the column of W which corresponds to column K of A */

#line 248 "clahef.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 252 "clahef.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 252 "clahef.f"
	    goto L30;
#line 252 "clahef.f"
	}

#line 255 "clahef.f"
	kstep = 1;

/*        Copy column K of A to column KW of W and update it */

#line 259 "clahef.f"
	i__1 = k - 1;
#line 259 "clahef.f"
	ccopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
#line 260 "clahef.f"
	i__1 = k + kw * w_dim1;
#line 260 "clahef.f"
	i__2 = k + k * a_dim1;
#line 260 "clahef.f"
	d__1 = a[i__2].r;
#line 260 "clahef.f"
	w[i__1].r = d__1, w[i__1].i = 0.;
#line 261 "clahef.f"
	if (k < *n) {
#line 262 "clahef.f"
	    i__1 = *n - k;
#line 262 "clahef.f"
	    z__1.r = -1., z__1.i = -0.;
#line 262 "clahef.f"
	    cgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 264 "clahef.f"
	    i__1 = k + kw * w_dim1;
#line 264 "clahef.f"
	    i__2 = k + kw * w_dim1;
#line 264 "clahef.f"
	    d__1 = w[i__2].r;
#line 264 "clahef.f"
	    w[i__1].r = d__1, w[i__1].i = 0.;
#line 265 "clahef.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 270 "clahef.f"
	i__1 = k + kw * w_dim1;
#line 270 "clahef.f"
	absakk = (d__1 = w[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 276 "clahef.f"
	if (k > 1) {
#line 277 "clahef.f"
	    i__1 = k - 1;
#line 277 "clahef.f"
	    imax = icamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 278 "clahef.f"
	    i__1 = imax + kw * w_dim1;
#line 278 "clahef.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    kw * w_dim1]), abs(d__2));
#line 279 "clahef.f"
	} else {
#line 280 "clahef.f"
	    colmax = 0.;
#line 281 "clahef.f"
	}

#line 283 "clahef.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 287 "clahef.f"
	    if (*info == 0) {
#line 287 "clahef.f"
		*info = k;
#line 287 "clahef.f"
	    }
#line 289 "clahef.f"
	    kp = k;
#line 290 "clahef.f"
	    i__1 = k + k * a_dim1;
#line 290 "clahef.f"
	    i__2 = k + k * a_dim1;
#line 290 "clahef.f"
	    d__1 = a[i__2].r;
#line 290 "clahef.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 291 "clahef.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
#line 298 "clahef.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 302 "clahef.f"
		kp = k;
#line 303 "clahef.f"
	    } else {

/*              BEGIN pivot search along IMAX row */


/*              Copy column IMAX to column KW-1 of W and update it */

#line 310 "clahef.f"
		i__1 = imax - 1;
#line 310 "clahef.f"
		ccopy_(&i__1, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
#line 311 "clahef.f"
		i__1 = imax + (kw - 1) * w_dim1;
#line 311 "clahef.f"
		i__2 = imax + imax * a_dim1;
#line 311 "clahef.f"
		d__1 = a[i__2].r;
#line 311 "clahef.f"
		w[i__1].r = d__1, w[i__1].i = 0.;
#line 312 "clahef.f"
		i__1 = k - imax;
#line 312 "clahef.f"
		ccopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);
#line 314 "clahef.f"
		i__1 = k - imax;
#line 314 "clahef.f"
		clacgv_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);
#line 315 "clahef.f"
		if (k < *n) {
#line 316 "clahef.f"
		    i__1 = *n - k;
#line 316 "clahef.f"
		    z__1.r = -1., z__1.i = -0.;
#line 316 "clahef.f"
		    cgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 319 "clahef.f"
		    i__1 = imax + (kw - 1) * w_dim1;
#line 319 "clahef.f"
		    i__2 = imax + (kw - 1) * w_dim1;
#line 319 "clahef.f"
		    d__1 = w[i__2].r;
#line 319 "clahef.f"
		    w[i__1].r = d__1, w[i__1].i = 0.;
#line 320 "clahef.f"
		}

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value. */
/*              Determine only ROWMAX. */

#line 326 "clahef.f"
		i__1 = k - imax;
#line 326 "clahef.f"
		jmax = imax + icamax_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1],
			 &c__1);
#line 327 "clahef.f"
		i__1 = jmax + (kw - 1) * w_dim1;
#line 327 "clahef.f"
		rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			jmax + (kw - 1) * w_dim1]), abs(d__2));
#line 328 "clahef.f"
		if (imax > 1) {
#line 329 "clahef.f"
		    i__1 = imax - 1;
#line 329 "clahef.f"
		    jmax = icamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
/* Computing MAX */
#line 330 "clahef.f"
		    i__1 = jmax + (kw - 1) * w_dim1;
#line 330 "clahef.f"
		    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&w[jmax + (kw - 1) * w_dim1]), abs(
			    d__2));
#line 330 "clahef.f"
		    rowmax = max(d__3,d__4);
#line 331 "clahef.f"
		}

/*              Case(2) */
#line 334 "clahef.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 338 "clahef.f"
		    kp = k;

/*              Case(3) */
#line 341 "clahef.f"
		} else /* if(complicated condition) */ {
#line 341 "clahef.f"
		    i__1 = imax + (kw - 1) * w_dim1;
#line 341 "clahef.f"
		    if ((d__1 = w[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 347 "clahef.f"
			kp = imax;

/*                 copy column KW-1 of W to column KW of W */

#line 351 "clahef.f"
			ccopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
				w_dim1 + 1], &c__1);

/*              Case(4) */
#line 354 "clahef.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 359 "clahef.f"
			kp = imax;
#line 360 "clahef.f"
			kstep = 2;
#line 361 "clahef.f"
		    }
#line 361 "clahef.f"
		}


/*              END pivot search along IMAX row */

#line 366 "clahef.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 374 "clahef.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 378 "clahef.f"
	    kkw = *nb + kk - *n;

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KKW of W. */

#line 383 "clahef.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K-1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 390 "clahef.f"
		i__1 = kp + kp * a_dim1;
#line 390 "clahef.f"
		i__2 = kk + kk * a_dim1;
#line 390 "clahef.f"
		d__1 = a[i__2].r;
#line 390 "clahef.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 391 "clahef.f"
		i__1 = kk - 1 - kp;
#line 391 "clahef.f"
		ccopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 393 "clahef.f"
		i__1 = kk - 1 - kp;
#line 393 "clahef.f"
		clacgv_(&i__1, &a[kp + (kp + 1) * a_dim1], lda);
#line 394 "clahef.f"
		if (kp > 1) {
#line 394 "clahef.f"
		    i__1 = kp - 1;
#line 394 "clahef.f"
		    ccopy_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 394 "clahef.f"
		}

/*              Interchange rows KK and KP in last K+1 to N columns of A */
/*              (columns K (or K and K-1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in last KKW to NB columns of W. */

#line 402 "clahef.f"
		if (k < *n) {
#line 402 "clahef.f"
		    i__1 = *n - k;
#line 402 "clahef.f"
		    cswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 402 "clahef.f"
		}
#line 405 "clahef.f"
		i__1 = *n - kk + 1;
#line 405 "clahef.f"
		cswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 407 "clahef.f"
	    }

#line 409 "clahef.f"
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
/*              A( K, K ) = DBLE( W( K, K) ) to separately copy diagonal */
/*              element D(k,k) from W (potentially saves only one load)) */
#line 427 "clahef.f"
		ccopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 428 "clahef.f"
		if (k > 1) {

/*                 (NOTE: No need to check if A(k,k) is NOT ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  case A(k,k) = 0 falls into 2x2 pivot case(4)) */

#line 434 "clahef.f"
		    i__1 = k + k * a_dim1;
#line 434 "clahef.f"
		    r1 = 1. / a[i__1].r;
#line 435 "clahef.f"
		    i__1 = k - 1;
#line 435 "clahef.f"
		    csscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);

/*                 (2) Conjugate column W(kw) */

#line 439 "clahef.f"
		    i__1 = k - 1;
#line 439 "clahef.f"
		    clacgv_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 440 "clahef.f"
		}

#line 442 "clahef.f"
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

#line 459 "clahef.f"
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

/*                 = ( conj(D21)*( D11 ) D21*(  -1 ) ) */
/*                   (           (  -1 )     ( D22 ) ), */

/*                 where D11 = d22/d21, */
/*                       D22 = d11/conj(d21), */
/*                       D21 = T/d21, */
/*                       T = 1/(D22*D11-1). */

/*                 (NOTE: No need to check for division by ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  (a) d21 != 0, since in 2x2 pivot case(4) */
/*                      |d21| should be larger than |d11| and |d22|; */
/*                  (b) (D22*D11 - 1) != 0, since from (a), */
/*                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.) */

#line 503 "clahef.f"
		    i__1 = k - 1 + kw * w_dim1;
#line 503 "clahef.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 504 "clahef.f"
		    d_cnjg(&z__2, &d21);
#line 504 "clahef.f"
		    z_div(&z__1, &w[k + kw * w_dim1], &z__2);
#line 504 "clahef.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 505 "clahef.f"
		    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d21);
#line 505 "clahef.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 506 "clahef.f"
		    z__1.r = d11.r * d22.r - d11.i * d22.i, z__1.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 506 "clahef.f"
		    t = 1. / (z__1.r - 1.);
#line 507 "clahef.f"
		    z__2.r = t, z__2.i = 0.;
#line 507 "clahef.f"
		    z_div(&z__1, &z__2, &d21);
#line 507 "clahef.f"
		    d21.r = z__1.r, d21.i = z__1.i;

/*                 Update elements in columns A(k-1) and A(k) as */
/*                 dot products of rows of ( W(kw-1) W(kw) ) and columns */
/*                 of D**(-1) */

#line 513 "clahef.f"
		    i__1 = k - 2;
#line 513 "clahef.f"
		    for (j = 1; j <= i__1; ++j) {
#line 514 "clahef.f"
			i__2 = j + (k - 1) * a_dim1;
#line 514 "clahef.f"
			i__3 = j + (kw - 1) * w_dim1;
#line 514 "clahef.f"
			z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__3.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 514 "clahef.f"
			i__4 = j + kw * w_dim1;
#line 514 "clahef.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 514 "clahef.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 514 "clahef.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 515 "clahef.f"
			i__2 = j + k * a_dim1;
#line 515 "clahef.f"
			d_cnjg(&z__2, &d21);
#line 515 "clahef.f"
			i__3 = j + kw * w_dim1;
#line 515 "clahef.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 515 "clahef.f"
			i__4 = j + (kw - 1) * w_dim1;
#line 515 "clahef.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 515 "clahef.f"
			z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
				z__2.r * z__3.i + z__2.i * z__3.r;
#line 515 "clahef.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 517 "clahef.f"
/* L20: */
#line 517 "clahef.f"
		    }
#line 518 "clahef.f"
		}

/*              Copy D(k) to A */

#line 522 "clahef.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 522 "clahef.f"
		i__2 = k - 1 + (kw - 1) * w_dim1;
#line 522 "clahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 523 "clahef.f"
		i__1 = k - 1 + k * a_dim1;
#line 523 "clahef.f"
		i__2 = k - 1 + kw * w_dim1;
#line 523 "clahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 524 "clahef.f"
		i__1 = k + k * a_dim1;
#line 524 "clahef.f"
		i__2 = k + kw * w_dim1;
#line 524 "clahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

/*              (2) Conjugate columns W(kw) and W(kw-1) */

#line 528 "clahef.f"
		i__1 = k - 1;
#line 528 "clahef.f"
		clacgv_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 529 "clahef.f"
		i__1 = k - 2;
#line 529 "clahef.f"
		clacgv_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);

#line 531 "clahef.f"
	    }

#line 533 "clahef.f"
	}

/*        Store details of the interchanges in IPIV */

#line 537 "clahef.f"
	if (kstep == 1) {
#line 538 "clahef.f"
	    ipiv[k] = kp;
#line 539 "clahef.f"
	} else {
#line 540 "clahef.f"
	    ipiv[k] = -kp;
#line 541 "clahef.f"
	    ipiv[k - 1] = -kp;
#line 542 "clahef.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 546 "clahef.f"
	k -= kstep;
#line 547 "clahef.f"
	goto L10;

#line 549 "clahef.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**H = A11 - U12*W**H */

/*        computing blocks of NB columns at a time (note that conjg(W) is */
/*        actually stored) */

#line 558 "clahef.f"
	i__1 = -(*nb);
#line 558 "clahef.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 559 "clahef.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 559 "clahef.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 563 "clahef.f"
	    i__2 = j + jb - 1;
#line 563 "clahef.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 564 "clahef.f"
		i__3 = jj + jj * a_dim1;
#line 564 "clahef.f"
		i__4 = jj + jj * a_dim1;
#line 564 "clahef.f"
		d__1 = a[i__4].r;
#line 564 "clahef.f"
		a[i__3].r = d__1, a[i__3].i = 0.;
#line 565 "clahef.f"
		i__3 = jj - j + 1;
#line 565 "clahef.f"
		i__4 = *n - k;
#line 565 "clahef.f"
		z__1.r = -1., z__1.i = -0.;
#line 565 "clahef.f"
		cgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, 
			&a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 568 "clahef.f"
		i__3 = jj + jj * a_dim1;
#line 568 "clahef.f"
		i__4 = jj + jj * a_dim1;
#line 568 "clahef.f"
		d__1 = a[i__4].r;
#line 568 "clahef.f"
		a[i__3].r = d__1, a[i__3].i = 0.;
#line 569 "clahef.f"
/* L40: */
#line 569 "clahef.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 573 "clahef.f"
	    i__2 = j - 1;
#line 573 "clahef.f"
	    i__3 = *n - k;
#line 573 "clahef.f"
	    z__1.r = -1., z__1.i = -0.;
#line 573 "clahef.f"
	    cgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, &a[(
		    k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * w_dim1], ldw,
		     &c_b1, &a[j * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
#line 576 "clahef.f"
/* L50: */
#line 576 "clahef.f"
	}

/*        Put U12 in standard form by partially undoing the interchanges */
/*        in of rows in columns k+1:n looping backwards from k+1 to n */

#line 581 "clahef.f"
	j = k + 1;
#line 582 "clahef.f"
L60:

/*           Undo the interchanges (if any) of rows J and JP */
/*           at each step J */

/*           (Here, J is a diagonal index) */
#line 588 "clahef.f"
	jj = j;
#line 589 "clahef.f"
	jp = ipiv[j];
#line 590 "clahef.f"
	if (jp < 0) {
#line 591 "clahef.f"
	    jp = -jp;
/*              (Here, J is a diagonal index) */
#line 593 "clahef.f"
	    ++j;
#line 594 "clahef.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length N-J+1 */
/*           of the rows to swap back doesn't include diagonal element) */
#line 597 "clahef.f"
	++j;
#line 598 "clahef.f"
	if (jp != jj && j <= *n) {
#line 598 "clahef.f"
	    i__1 = *n - j + 1;
#line 598 "clahef.f"
	    cswap_(&i__1, &a[jp + j * a_dim1], lda, &a[jj + j * a_dim1], lda);
#line 598 "clahef.f"
	}
#line 600 "clahef.f"
	if (j <= *n) {
#line 600 "clahef.f"
	    goto L60;
#line 600 "clahef.f"
	}

/*        Set KB to the number of columns factorized */

#line 605 "clahef.f"
	*kb = *n - k;

#line 607 "clahef.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 (note that conjg(W) is actually stored) */

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 615 "clahef.f"
	k = 1;
#line 616 "clahef.f"
L70:

/*        Exit from loop */

#line 620 "clahef.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 620 "clahef.f"
	    goto L90;
#line 620 "clahef.f"
	}

#line 623 "clahef.f"
	kstep = 1;

/*        Copy column K of A to column K of W and update it */

#line 627 "clahef.f"
	i__1 = k + k * w_dim1;
#line 627 "clahef.f"
	i__2 = k + k * a_dim1;
#line 627 "clahef.f"
	d__1 = a[i__2].r;
#line 627 "clahef.f"
	w[i__1].r = d__1, w[i__1].i = 0.;
#line 628 "clahef.f"
	if (k < *n) {
#line 628 "clahef.f"
	    i__1 = *n - k;
#line 628 "clahef.f"
	    ccopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &w[k + 1 + k * 
		    w_dim1], &c__1);
#line 628 "clahef.f"
	}
#line 630 "clahef.f"
	i__1 = *n - k + 1;
#line 630 "clahef.f"
	i__2 = k - 1;
#line 630 "clahef.f"
	z__1.r = -1., z__1.i = -0.;
#line 630 "clahef.f"
	cgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[k 
		+ w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1, (ftnlen)12);
#line 632 "clahef.f"
	i__1 = k + k * w_dim1;
#line 632 "clahef.f"
	i__2 = k + k * w_dim1;
#line 632 "clahef.f"
	d__1 = w[i__2].r;
#line 632 "clahef.f"
	w[i__1].r = d__1, w[i__1].i = 0.;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 637 "clahef.f"
	i__1 = k + k * w_dim1;
#line 637 "clahef.f"
	absakk = (d__1 = w[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 643 "clahef.f"
	if (k < *n) {
#line 644 "clahef.f"
	    i__1 = *n - k;
#line 644 "clahef.f"
	    imax = k + icamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 645 "clahef.f"
	    i__1 = imax + k * w_dim1;
#line 645 "clahef.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    k * w_dim1]), abs(d__2));
#line 646 "clahef.f"
	} else {
#line 647 "clahef.f"
	    colmax = 0.;
#line 648 "clahef.f"
	}

#line 650 "clahef.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 654 "clahef.f"
	    if (*info == 0) {
#line 654 "clahef.f"
		*info = k;
#line 654 "clahef.f"
	    }
#line 656 "clahef.f"
	    kp = k;
#line 657 "clahef.f"
	    i__1 = k + k * a_dim1;
#line 657 "clahef.f"
	    i__2 = k + k * a_dim1;
#line 657 "clahef.f"
	    d__1 = a[i__2].r;
#line 657 "clahef.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 658 "clahef.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
#line 665 "clahef.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 669 "clahef.f"
		kp = k;
#line 670 "clahef.f"
	    } else {

/*              BEGIN pivot search along IMAX row */


/*              Copy column IMAX to column K+1 of W and update it */

#line 677 "clahef.f"
		i__1 = imax - k;
#line 677 "clahef.f"
		ccopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 678 "clahef.f"
		i__1 = imax - k;
#line 678 "clahef.f"
		clacgv_(&i__1, &w[k + (k + 1) * w_dim1], &c__1);
#line 679 "clahef.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 679 "clahef.f"
		i__2 = imax + imax * a_dim1;
#line 679 "clahef.f"
		d__1 = a[i__2].r;
#line 679 "clahef.f"
		w[i__1].r = d__1, w[i__1].i = 0.;
#line 680 "clahef.f"
		if (imax < *n) {
#line 680 "clahef.f"
		    i__1 = *n - imax;
#line 680 "clahef.f"
		    ccopy_(&i__1, &a[imax + 1 + imax * a_dim1], &c__1, &w[
			    imax + 1 + (k + 1) * w_dim1], &c__1);
#line 680 "clahef.f"
		}
#line 683 "clahef.f"
		i__1 = *n - k + 1;
#line 683 "clahef.f"
		i__2 = k - 1;
#line 683 "clahef.f"
		z__1.r = -1., z__1.i = -0.;
#line 683 "clahef.f"
		cgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], 
			lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 1) * 
			w_dim1], &c__1, (ftnlen)12);
#line 686 "clahef.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 686 "clahef.f"
		i__2 = imax + (k + 1) * w_dim1;
#line 686 "clahef.f"
		d__1 = w[i__2].r;
#line 686 "clahef.f"
		w[i__1].r = d__1, w[i__1].i = 0.;

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value. */
/*              Determine only ROWMAX. */

#line 692 "clahef.f"
		i__1 = imax - k;
#line 692 "clahef.f"
		jmax = k - 1 + icamax_(&i__1, &w[k + (k + 1) * w_dim1], &c__1)
			;
#line 693 "clahef.f"
		i__1 = jmax + (k + 1) * w_dim1;
#line 693 "clahef.f"
		rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			jmax + (k + 1) * w_dim1]), abs(d__2));
#line 694 "clahef.f"
		if (imax < *n) {
#line 695 "clahef.f"
		    i__1 = *n - imax;
#line 695 "clahef.f"
		    jmax = imax + icamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
/* Computing MAX */
#line 696 "clahef.f"
		    i__1 = jmax + (k + 1) * w_dim1;
#line 696 "clahef.f"
		    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&w[jmax + (k + 1) * w_dim1]), abs(
			    d__2));
#line 696 "clahef.f"
		    rowmax = max(d__3,d__4);
#line 697 "clahef.f"
		}

/*              Case(2) */
#line 700 "clahef.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 704 "clahef.f"
		    kp = k;

/*              Case(3) */
#line 707 "clahef.f"
		} else /* if(complicated condition) */ {
#line 707 "clahef.f"
		    i__1 = imax + (k + 1) * w_dim1;
#line 707 "clahef.f"
		    if ((d__1 = w[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 713 "clahef.f"
			kp = imax;

/*                 copy column K+1 of W to column K of W */

#line 717 "clahef.f"
			i__1 = *n - k + 1;
#line 717 "clahef.f"
			ccopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + 
				k * w_dim1], &c__1);

/*              Case(4) */
#line 720 "clahef.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 725 "clahef.f"
			kp = imax;
#line 726 "clahef.f"
			kstep = 2;
#line 727 "clahef.f"
		    }
#line 727 "clahef.f"
		}


/*              END pivot search along IMAX row */

#line 732 "clahef.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 740 "clahef.f"
	    kk = k + kstep - 1;

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KK of W. */

#line 745 "clahef.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K+1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 752 "clahef.f"
		i__1 = kp + kp * a_dim1;
#line 752 "clahef.f"
		i__2 = kk + kk * a_dim1;
#line 752 "clahef.f"
		d__1 = a[i__2].r;
#line 752 "clahef.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 753 "clahef.f"
		i__1 = kp - kk - 1;
#line 753 "clahef.f"
		ccopy_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 755 "clahef.f"
		i__1 = kp - kk - 1;
#line 755 "clahef.f"
		clacgv_(&i__1, &a[kp + (kk + 1) * a_dim1], lda);
#line 756 "clahef.f"
		if (kp < *n) {
#line 756 "clahef.f"
		    i__1 = *n - kp;
#line 756 "clahef.f"
		    ccopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 756 "clahef.f"
		}

/*              Interchange rows KK and KP in first K-1 columns of A */
/*              (columns K (or K and K+1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in first KK columns of W. */

#line 764 "clahef.f"
		if (k > 1) {
#line 764 "clahef.f"
		    i__1 = k - 1;
#line 764 "clahef.f"
		    cswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 764 "clahef.f"
		}
#line 766 "clahef.f"
		cswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 767 "clahef.f"
	    }

#line 769 "clahef.f"
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
/*              A( K, K ) = DBLE( W( K, K) ) to separately copy diagonal */
/*              element D(k,k) from W (potentially saves only one load)) */
#line 787 "clahef.f"
		i__1 = *n - k + 1;
#line 787 "clahef.f"
		ccopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 788 "clahef.f"
		if (k < *n) {

/*                 (NOTE: No need to check if A(k,k) is NOT ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  case A(k,k) = 0 falls into 2x2 pivot case(4)) */

#line 794 "clahef.f"
		    i__1 = k + k * a_dim1;
#line 794 "clahef.f"
		    r1 = 1. / a[i__1].r;
#line 795 "clahef.f"
		    i__1 = *n - k;
#line 795 "clahef.f"
		    csscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);

/*                 (2) Conjugate column W(k) */

#line 799 "clahef.f"
		    i__1 = *n - k;
#line 799 "clahef.f"
		    clacgv_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 800 "clahef.f"
		}

#line 802 "clahef.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

/*              (1) Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2 */
/*              block D(k:k+1,k:k+1) in columns k and k+1 of A. */
/*              (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT */
/*              block and not stored) */
/*                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1) */
/*                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) = */
/*                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) ) */

#line 819 "clahef.f"
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

/*                 = ( conj(D21)*( D11 ) D21*(  -1 ) ) */
/*                   (           (  -1 )     ( D22 ) ) */

/*                 where D11 = d22/d21, */
/*                       D22 = d11/conj(d21), */
/*                       D21 = T/d21, */
/*                       T = 1/(D22*D11-1). */

/*                 (NOTE: No need to check for division by ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  (a) d21 != 0, since in 2x2 pivot case(4) */
/*                      |d21| should be larger than |d11| and |d22|; */
/*                  (b) (D22*D11 - 1) != 0, since from (a), */
/*                      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.) */

#line 863 "clahef.f"
		    i__1 = k + 1 + k * w_dim1;
#line 863 "clahef.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 864 "clahef.f"
		    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
#line 864 "clahef.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 865 "clahef.f"
		    d_cnjg(&z__2, &d21);
#line 865 "clahef.f"
		    z_div(&z__1, &w[k + k * w_dim1], &z__2);
#line 865 "clahef.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 866 "clahef.f"
		    z__1.r = d11.r * d22.r - d11.i * d22.i, z__1.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 866 "clahef.f"
		    t = 1. / (z__1.r - 1.);
#line 867 "clahef.f"
		    z__2.r = t, z__2.i = 0.;
#line 867 "clahef.f"
		    z_div(&z__1, &z__2, &d21);
#line 867 "clahef.f"
		    d21.r = z__1.r, d21.i = z__1.i;

/*                 Update elements in columns A(k) and A(k+1) as */
/*                 dot products of rows of ( W(k) W(k+1) ) and columns */
/*                 of D**(-1) */

#line 873 "clahef.f"
		    i__1 = *n;
#line 873 "clahef.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 874 "clahef.f"
			i__2 = j + k * a_dim1;
#line 874 "clahef.f"
			d_cnjg(&z__2, &d21);
#line 874 "clahef.f"
			i__3 = j + k * w_dim1;
#line 874 "clahef.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 874 "clahef.f"
			i__4 = j + (k + 1) * w_dim1;
#line 874 "clahef.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 874 "clahef.f"
			z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
				z__2.r * z__3.i + z__2.i * z__3.r;
#line 874 "clahef.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 876 "clahef.f"
			i__2 = j + (k + 1) * a_dim1;
#line 876 "clahef.f"
			i__3 = j + (k + 1) * w_dim1;
#line 876 "clahef.f"
			z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__3.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 876 "clahef.f"
			i__4 = j + k * w_dim1;
#line 876 "clahef.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 876 "clahef.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 876 "clahef.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 877 "clahef.f"
/* L80: */
#line 877 "clahef.f"
		    }
#line 878 "clahef.f"
		}

/*              Copy D(k) to A */

#line 882 "clahef.f"
		i__1 = k + k * a_dim1;
#line 882 "clahef.f"
		i__2 = k + k * w_dim1;
#line 882 "clahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 883 "clahef.f"
		i__1 = k + 1 + k * a_dim1;
#line 883 "clahef.f"
		i__2 = k + 1 + k * w_dim1;
#line 883 "clahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 884 "clahef.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 884 "clahef.f"
		i__2 = k + 1 + (k + 1) * w_dim1;
#line 884 "clahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

/*              (2) Conjugate columns W(k) and W(k+1) */

#line 888 "clahef.f"
		i__1 = *n - k;
#line 888 "clahef.f"
		clacgv_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 889 "clahef.f"
		i__1 = *n - k - 1;
#line 889 "clahef.f"
		clacgv_(&i__1, &w[k + 2 + (k + 1) * w_dim1], &c__1);

#line 891 "clahef.f"
	    }

#line 893 "clahef.f"
	}

/*        Store details of the interchanges in IPIV */

#line 897 "clahef.f"
	if (kstep == 1) {
#line 898 "clahef.f"
	    ipiv[k] = kp;
#line 899 "clahef.f"
	} else {
#line 900 "clahef.f"
	    ipiv[k] = -kp;
#line 901 "clahef.f"
	    ipiv[k + 1] = -kp;
#line 902 "clahef.f"
	}

/*        Increase K and return to the start of the main loop */

#line 906 "clahef.f"
	k += kstep;
#line 907 "clahef.f"
	goto L70;

#line 909 "clahef.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**H = A22 - L21*W**H */

/*        computing blocks of NB columns at a time (note that conjg(W) is */
/*        actually stored) */

#line 918 "clahef.f"
	i__1 = *n;
#line 918 "clahef.f"
	i__2 = *nb;
#line 918 "clahef.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 919 "clahef.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 919 "clahef.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 923 "clahef.f"
	    i__3 = j + jb - 1;
#line 923 "clahef.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 924 "clahef.f"
		i__4 = jj + jj * a_dim1;
#line 924 "clahef.f"
		i__5 = jj + jj * a_dim1;
#line 924 "clahef.f"
		d__1 = a[i__5].r;
#line 924 "clahef.f"
		a[i__4].r = d__1, a[i__4].i = 0.;
#line 925 "clahef.f"
		i__4 = j + jb - jj;
#line 925 "clahef.f"
		i__5 = k - 1;
#line 925 "clahef.f"
		z__1.r = -1., z__1.i = -0.;
#line 925 "clahef.f"
		cgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1]
			, &c__1, (ftnlen)12);
#line 928 "clahef.f"
		i__4 = jj + jj * a_dim1;
#line 928 "clahef.f"
		i__5 = jj + jj * a_dim1;
#line 928 "clahef.f"
		d__1 = a[i__5].r;
#line 928 "clahef.f"
		a[i__4].r = d__1, a[i__4].i = 0.;
#line 929 "clahef.f"
/* L100: */
#line 929 "clahef.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 933 "clahef.f"
	    if (j + jb <= *n) {
#line 933 "clahef.f"
		i__3 = *n - j - jb + 1;
#line 933 "clahef.f"
		i__4 = k - 1;
#line 933 "clahef.f"
		z__1.r = -1., z__1.i = -0.;
#line 933 "clahef.f"
		cgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, 
			&a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 933 "clahef.f"
	    }
#line 937 "clahef.f"
/* L110: */
#line 937 "clahef.f"
	}

/*        Put L21 in standard form by partially undoing the interchanges */
/*        of rows in columns 1:k-1 looping backwards from k-1 to 1 */

#line 942 "clahef.f"
	j = k - 1;
#line 943 "clahef.f"
L120:

/*           Undo the interchanges (if any) of rows J and JP */
/*           at each step J */

/*           (Here, J is a diagonal index) */
#line 949 "clahef.f"
	jj = j;
#line 950 "clahef.f"
	jp = ipiv[j];
#line 951 "clahef.f"
	if (jp < 0) {
#line 952 "clahef.f"
	    jp = -jp;
/*              (Here, J is a diagonal index) */
#line 954 "clahef.f"
	    --j;
#line 955 "clahef.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length J */
/*           of the rows to swap back doesn't include diagonal element) */
#line 958 "clahef.f"
	--j;
#line 959 "clahef.f"
	if (jp != jj && j >= 1) {
#line 959 "clahef.f"
	    cswap_(&j, &a[jp + a_dim1], lda, &a[jj + a_dim1], lda);
#line 959 "clahef.f"
	}
#line 961 "clahef.f"
	if (j >= 1) {
#line 961 "clahef.f"
	    goto L120;
#line 961 "clahef.f"
	}

/*        Set KB to the number of columns factorized */

#line 966 "clahef.f"
	*kb = k - 1;

#line 968 "clahef.f"
    }
#line 969 "clahef.f"
    return 0;

/*     End of CLAHEF */

} /* clahef_ */

