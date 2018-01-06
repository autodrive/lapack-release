#line 1 "zlahef.f"
/* zlahef.f -- translated by f2c (version 20100827).
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

#line 1 "zlahef.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLAHEF computes a partial factorization of a complex Hermitian indefinite matrix using the Bunc
h-Kaufman diagonal pivoting method (blocked algorithm, calling Level 3 BLAS). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAHEF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahef.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahef.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahef.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAHEF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */

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
/* > ZLAHEF computes a partial factorization of a complex Hermitian */
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
/* > ZLAHEF is an auxiliary routine called by ZHETRF. It uses blocked code */
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

/* > \date December 2016 */

/* > \ingroup complex16HEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  December 2016,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zlahef_(char *uplo, integer *n, integer *nb, integer *kb,
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
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
    static doublereal absakk;
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal colmax;
    extern /* Subroutine */ int zlacgv_(integer *, doublecomplex *, integer *)
	    ;
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

#line 229 "zlahef.f"
    /* Parameter adjustments */
#line 229 "zlahef.f"
    a_dim1 = *lda;
#line 229 "zlahef.f"
    a_offset = 1 + a_dim1;
#line 229 "zlahef.f"
    a -= a_offset;
#line 229 "zlahef.f"
    --ipiv;
#line 229 "zlahef.f"
    w_dim1 = *ldw;
#line 229 "zlahef.f"
    w_offset = 1 + w_dim1;
#line 229 "zlahef.f"
    w -= w_offset;
#line 229 "zlahef.f"

#line 229 "zlahef.f"
    /* Function Body */
#line 229 "zlahef.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 233 "zlahef.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 235 "zlahef.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 (note that conjg(W) is actually stored) */

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

/*        KW is the column of W which corresponds to column K of A */

#line 245 "zlahef.f"
	k = *n;
#line 246 "zlahef.f"
L10:
#line 247 "zlahef.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 251 "zlahef.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 251 "zlahef.f"
	    goto L30;
#line 251 "zlahef.f"
	}

#line 254 "zlahef.f"
	kstep = 1;

/*        Copy column K of A to column KW of W and update it */

#line 258 "zlahef.f"
	i__1 = k - 1;
#line 258 "zlahef.f"
	zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
#line 259 "zlahef.f"
	i__1 = k + kw * w_dim1;
#line 259 "zlahef.f"
	i__2 = k + k * a_dim1;
#line 259 "zlahef.f"
	d__1 = a[i__2].r;
#line 259 "zlahef.f"
	w[i__1].r = d__1, w[i__1].i = 0.;
#line 260 "zlahef.f"
	if (k < *n) {
#line 261 "zlahef.f"
	    i__1 = *n - k;
#line 261 "zlahef.f"
	    z__1.r = -1., z__1.i = -0.;
#line 261 "zlahef.f"
	    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 263 "zlahef.f"
	    i__1 = k + kw * w_dim1;
#line 263 "zlahef.f"
	    i__2 = k + kw * w_dim1;
#line 263 "zlahef.f"
	    d__1 = w[i__2].r;
#line 263 "zlahef.f"
	    w[i__1].r = d__1, w[i__1].i = 0.;
#line 264 "zlahef.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 269 "zlahef.f"
	i__1 = k + kw * w_dim1;
#line 269 "zlahef.f"
	absakk = (d__1 = w[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 275 "zlahef.f"
	if (k > 1) {
#line 276 "zlahef.f"
	    i__1 = k - 1;
#line 276 "zlahef.f"
	    imax = izamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 277 "zlahef.f"
	    i__1 = imax + kw * w_dim1;
#line 277 "zlahef.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    kw * w_dim1]), abs(d__2));
#line 278 "zlahef.f"
	} else {
#line 279 "zlahef.f"
	    colmax = 0.;
#line 280 "zlahef.f"
	}

#line 282 "zlahef.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 286 "zlahef.f"
	    if (*info == 0) {
#line 286 "zlahef.f"
		*info = k;
#line 286 "zlahef.f"
	    }
#line 288 "zlahef.f"
	    kp = k;
#line 289 "zlahef.f"
	    i__1 = k + k * a_dim1;
#line 289 "zlahef.f"
	    i__2 = k + k * a_dim1;
#line 289 "zlahef.f"
	    d__1 = a[i__2].r;
#line 289 "zlahef.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 290 "zlahef.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
#line 297 "zlahef.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 301 "zlahef.f"
		kp = k;
#line 302 "zlahef.f"
	    } else {

/*              BEGIN pivot search along IMAX row */


/*              Copy column IMAX to column KW-1 of W and update it */

#line 309 "zlahef.f"
		i__1 = imax - 1;
#line 309 "zlahef.f"
		zcopy_(&i__1, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
#line 310 "zlahef.f"
		i__1 = imax + (kw - 1) * w_dim1;
#line 310 "zlahef.f"
		i__2 = imax + imax * a_dim1;
#line 310 "zlahef.f"
		d__1 = a[i__2].r;
#line 310 "zlahef.f"
		w[i__1].r = d__1, w[i__1].i = 0.;
#line 311 "zlahef.f"
		i__1 = k - imax;
#line 311 "zlahef.f"
		zcopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);
#line 313 "zlahef.f"
		i__1 = k - imax;
#line 313 "zlahef.f"
		zlacgv_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);
#line 314 "zlahef.f"
		if (k < *n) {
#line 315 "zlahef.f"
		    i__1 = *n - k;
#line 315 "zlahef.f"
		    z__1.r = -1., z__1.i = -0.;
#line 315 "zlahef.f"
		    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 318 "zlahef.f"
		    i__1 = imax + (kw - 1) * w_dim1;
#line 318 "zlahef.f"
		    i__2 = imax + (kw - 1) * w_dim1;
#line 318 "zlahef.f"
		    d__1 = w[i__2].r;
#line 318 "zlahef.f"
		    w[i__1].r = d__1, w[i__1].i = 0.;
#line 319 "zlahef.f"
		}

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value. */
/*              Determine only ROWMAX. */

#line 325 "zlahef.f"
		i__1 = k - imax;
#line 325 "zlahef.f"
		jmax = imax + izamax_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1],
			 &c__1);
#line 326 "zlahef.f"
		i__1 = jmax + (kw - 1) * w_dim1;
#line 326 "zlahef.f"
		rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			jmax + (kw - 1) * w_dim1]), abs(d__2));
#line 327 "zlahef.f"
		if (imax > 1) {
#line 328 "zlahef.f"
		    i__1 = imax - 1;
#line 328 "zlahef.f"
		    jmax = izamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
/* Computing MAX */
#line 329 "zlahef.f"
		    i__1 = jmax + (kw - 1) * w_dim1;
#line 329 "zlahef.f"
		    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&w[jmax + (kw - 1) * w_dim1]), abs(
			    d__2));
#line 329 "zlahef.f"
		    rowmax = max(d__3,d__4);
#line 330 "zlahef.f"
		}

/*              Case(2) */
#line 333 "zlahef.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 337 "zlahef.f"
		    kp = k;

/*              Case(3) */
#line 340 "zlahef.f"
		} else /* if(complicated condition) */ {
#line 340 "zlahef.f"
		    i__1 = imax + (kw - 1) * w_dim1;
#line 340 "zlahef.f"
		    if ((d__1 = w[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 346 "zlahef.f"
			kp = imax;

/*                 copy column KW-1 of W to column KW of W */

#line 350 "zlahef.f"
			zcopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
				w_dim1 + 1], &c__1);

/*              Case(4) */
#line 353 "zlahef.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 358 "zlahef.f"
			kp = imax;
#line 359 "zlahef.f"
			kstep = 2;
#line 360 "zlahef.f"
		    }
#line 360 "zlahef.f"
		}


/*              END pivot search along IMAX row */

#line 365 "zlahef.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 373 "zlahef.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 377 "zlahef.f"
	    kkw = *nb + kk - *n;

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KKW of W. */

#line 382 "zlahef.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K-1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 389 "zlahef.f"
		i__1 = kp + kp * a_dim1;
#line 389 "zlahef.f"
		i__2 = kk + kk * a_dim1;
#line 389 "zlahef.f"
		d__1 = a[i__2].r;
#line 389 "zlahef.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 390 "zlahef.f"
		i__1 = kk - 1 - kp;
#line 390 "zlahef.f"
		zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 392 "zlahef.f"
		i__1 = kk - 1 - kp;
#line 392 "zlahef.f"
		zlacgv_(&i__1, &a[kp + (kp + 1) * a_dim1], lda);
#line 393 "zlahef.f"
		if (kp > 1) {
#line 393 "zlahef.f"
		    i__1 = kp - 1;
#line 393 "zlahef.f"
		    zcopy_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 393 "zlahef.f"
		}

/*              Interchange rows KK and KP in last K+1 to N columns of A */
/*              (columns K (or K and K-1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in last KKW to NB columns of W. */

#line 401 "zlahef.f"
		if (k < *n) {
#line 401 "zlahef.f"
		    i__1 = *n - k;
#line 401 "zlahef.f"
		    zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 401 "zlahef.f"
		}
#line 404 "zlahef.f"
		i__1 = *n - kk + 1;
#line 404 "zlahef.f"
		zswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 406 "zlahef.f"
	    }

#line 408 "zlahef.f"
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
#line 426 "zlahef.f"
		zcopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 427 "zlahef.f"
		if (k > 1) {

/*                 (NOTE: No need to check if A(k,k) is NOT ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  case A(k,k) = 0 falls into 2x2 pivot case(4)) */

#line 433 "zlahef.f"
		    i__1 = k + k * a_dim1;
#line 433 "zlahef.f"
		    r1 = 1. / a[i__1].r;
#line 434 "zlahef.f"
		    i__1 = k - 1;
#line 434 "zlahef.f"
		    zdscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);

/*                 (2) Conjugate column W(kw) */

#line 438 "zlahef.f"
		    i__1 = k - 1;
#line 438 "zlahef.f"
		    zlacgv_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 439 "zlahef.f"
		}

#line 441 "zlahef.f"
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

#line 458 "zlahef.f"
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

#line 502 "zlahef.f"
		    i__1 = k - 1 + kw * w_dim1;
#line 502 "zlahef.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 503 "zlahef.f"
		    d_cnjg(&z__2, &d21);
#line 503 "zlahef.f"
		    z_div(&z__1, &w[k + kw * w_dim1], &z__2);
#line 503 "zlahef.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 504 "zlahef.f"
		    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d21);
#line 504 "zlahef.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 505 "zlahef.f"
		    z__1.r = d11.r * d22.r - d11.i * d22.i, z__1.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 505 "zlahef.f"
		    t = 1. / (z__1.r - 1.);
#line 506 "zlahef.f"
		    z__2.r = t, z__2.i = 0.;
#line 506 "zlahef.f"
		    z_div(&z__1, &z__2, &d21);
#line 506 "zlahef.f"
		    d21.r = z__1.r, d21.i = z__1.i;

/*                 Update elements in columns A(k-1) and A(k) as */
/*                 dot products of rows of ( W(kw-1) W(kw) ) and columns */
/*                 of D**(-1) */

#line 512 "zlahef.f"
		    i__1 = k - 2;
#line 512 "zlahef.f"
		    for (j = 1; j <= i__1; ++j) {
#line 513 "zlahef.f"
			i__2 = j + (k - 1) * a_dim1;
#line 513 "zlahef.f"
			i__3 = j + (kw - 1) * w_dim1;
#line 513 "zlahef.f"
			z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__3.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 513 "zlahef.f"
			i__4 = j + kw * w_dim1;
#line 513 "zlahef.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 513 "zlahef.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 513 "zlahef.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 514 "zlahef.f"
			i__2 = j + k * a_dim1;
#line 514 "zlahef.f"
			d_cnjg(&z__2, &d21);
#line 514 "zlahef.f"
			i__3 = j + kw * w_dim1;
#line 514 "zlahef.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 514 "zlahef.f"
			i__4 = j + (kw - 1) * w_dim1;
#line 514 "zlahef.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 514 "zlahef.f"
			z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
				z__2.r * z__3.i + z__2.i * z__3.r;
#line 514 "zlahef.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 516 "zlahef.f"
/* L20: */
#line 516 "zlahef.f"
		    }
#line 517 "zlahef.f"
		}

/*              Copy D(k) to A */

#line 521 "zlahef.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 521 "zlahef.f"
		i__2 = k - 1 + (kw - 1) * w_dim1;
#line 521 "zlahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 522 "zlahef.f"
		i__1 = k - 1 + k * a_dim1;
#line 522 "zlahef.f"
		i__2 = k - 1 + kw * w_dim1;
#line 522 "zlahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 523 "zlahef.f"
		i__1 = k + k * a_dim1;
#line 523 "zlahef.f"
		i__2 = k + kw * w_dim1;
#line 523 "zlahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

/*              (2) Conjugate columns W(kw) and W(kw-1) */

#line 527 "zlahef.f"
		i__1 = k - 1;
#line 527 "zlahef.f"
		zlacgv_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 528 "zlahef.f"
		i__1 = k - 2;
#line 528 "zlahef.f"
		zlacgv_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);

#line 530 "zlahef.f"
	    }

#line 532 "zlahef.f"
	}

/*        Store details of the interchanges in IPIV */

#line 536 "zlahef.f"
	if (kstep == 1) {
#line 537 "zlahef.f"
	    ipiv[k] = kp;
#line 538 "zlahef.f"
	} else {
#line 539 "zlahef.f"
	    ipiv[k] = -kp;
#line 540 "zlahef.f"
	    ipiv[k - 1] = -kp;
#line 541 "zlahef.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 545 "zlahef.f"
	k -= kstep;
#line 546 "zlahef.f"
	goto L10;

#line 548 "zlahef.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**H = A11 - U12*W**H */

/*        computing blocks of NB columns at a time (note that conjg(W) is */
/*        actually stored) */

#line 557 "zlahef.f"
	i__1 = -(*nb);
#line 557 "zlahef.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 558 "zlahef.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 558 "zlahef.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 562 "zlahef.f"
	    i__2 = j + jb - 1;
#line 562 "zlahef.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 563 "zlahef.f"
		i__3 = jj + jj * a_dim1;
#line 563 "zlahef.f"
		i__4 = jj + jj * a_dim1;
#line 563 "zlahef.f"
		d__1 = a[i__4].r;
#line 563 "zlahef.f"
		a[i__3].r = d__1, a[i__3].i = 0.;
#line 564 "zlahef.f"
		i__3 = jj - j + 1;
#line 564 "zlahef.f"
		i__4 = *n - k;
#line 564 "zlahef.f"
		z__1.r = -1., z__1.i = -0.;
#line 564 "zlahef.f"
		zgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, 
			&a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 567 "zlahef.f"
		i__3 = jj + jj * a_dim1;
#line 567 "zlahef.f"
		i__4 = jj + jj * a_dim1;
#line 567 "zlahef.f"
		d__1 = a[i__4].r;
#line 567 "zlahef.f"
		a[i__3].r = d__1, a[i__3].i = 0.;
#line 568 "zlahef.f"
/* L40: */
#line 568 "zlahef.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 572 "zlahef.f"
	    i__2 = j - 1;
#line 572 "zlahef.f"
	    i__3 = *n - k;
#line 572 "zlahef.f"
	    z__1.r = -1., z__1.i = -0.;
#line 572 "zlahef.f"
	    zgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, &a[(
		    k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * w_dim1], ldw,
		     &c_b1, &a[j * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
#line 575 "zlahef.f"
/* L50: */
#line 575 "zlahef.f"
	}

/*        Put U12 in standard form by partially undoing the interchanges */
/*        in columns k+1:n looping backwards from k+1 to n */

#line 580 "zlahef.f"
	j = k + 1;
#line 581 "zlahef.f"
L60:

/*           Undo the interchanges (if any) of rows JJ and JP at each */
/*           step J */

/*           (Here, J is a diagonal index) */
#line 587 "zlahef.f"
	jj = j;
#line 588 "zlahef.f"
	jp = ipiv[j];
#line 589 "zlahef.f"
	if (jp < 0) {
#line 590 "zlahef.f"
	    jp = -jp;
/*              (Here, J is a diagonal index) */
#line 592 "zlahef.f"
	    ++j;
#line 593 "zlahef.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length N-J+1 */
/*           of the rows to swap back doesn't include diagonal element) */
#line 596 "zlahef.f"
	++j;
#line 597 "zlahef.f"
	if (jp != jj && j <= *n) {
#line 597 "zlahef.f"
	    i__1 = *n - j + 1;
#line 597 "zlahef.f"
	    zswap_(&i__1, &a[jp + j * a_dim1], lda, &a[jj + j * a_dim1], lda);
#line 597 "zlahef.f"
	}
#line 599 "zlahef.f"
	if (j < *n) {
#line 599 "zlahef.f"
	    goto L60;
#line 599 "zlahef.f"
	}

/*        Set KB to the number of columns factorized */

#line 604 "zlahef.f"
	*kb = *n - k;

#line 606 "zlahef.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 (note that conjg(W) is actually stored) */

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 614 "zlahef.f"
	k = 1;
#line 615 "zlahef.f"
L70:

/*        Exit from loop */

#line 619 "zlahef.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 619 "zlahef.f"
	    goto L90;
#line 619 "zlahef.f"
	}

#line 622 "zlahef.f"
	kstep = 1;

/*        Copy column K of A to column K of W and update it */

#line 626 "zlahef.f"
	i__1 = k + k * w_dim1;
#line 626 "zlahef.f"
	i__2 = k + k * a_dim1;
#line 626 "zlahef.f"
	d__1 = a[i__2].r;
#line 626 "zlahef.f"
	w[i__1].r = d__1, w[i__1].i = 0.;
#line 627 "zlahef.f"
	if (k < *n) {
#line 627 "zlahef.f"
	    i__1 = *n - k;
#line 627 "zlahef.f"
	    zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &w[k + 1 + k * 
		    w_dim1], &c__1);
#line 627 "zlahef.f"
	}
#line 629 "zlahef.f"
	i__1 = *n - k + 1;
#line 629 "zlahef.f"
	i__2 = k - 1;
#line 629 "zlahef.f"
	z__1.r = -1., z__1.i = -0.;
#line 629 "zlahef.f"
	zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[k 
		+ w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1, (ftnlen)12);
#line 631 "zlahef.f"
	i__1 = k + k * w_dim1;
#line 631 "zlahef.f"
	i__2 = k + k * w_dim1;
#line 631 "zlahef.f"
	d__1 = w[i__2].r;
#line 631 "zlahef.f"
	w[i__1].r = d__1, w[i__1].i = 0.;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 636 "zlahef.f"
	i__1 = k + k * w_dim1;
#line 636 "zlahef.f"
	absakk = (d__1 = w[i__1].r, abs(d__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 642 "zlahef.f"
	if (k < *n) {
#line 643 "zlahef.f"
	    i__1 = *n - k;
#line 643 "zlahef.f"
	    imax = k + izamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 644 "zlahef.f"
	    i__1 = imax + k * w_dim1;
#line 644 "zlahef.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    k * w_dim1]), abs(d__2));
#line 645 "zlahef.f"
	} else {
#line 646 "zlahef.f"
	    colmax = 0.;
#line 647 "zlahef.f"
	}

#line 649 "zlahef.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 653 "zlahef.f"
	    if (*info == 0) {
#line 653 "zlahef.f"
		*info = k;
#line 653 "zlahef.f"
	    }
#line 655 "zlahef.f"
	    kp = k;
#line 656 "zlahef.f"
	    i__1 = k + k * a_dim1;
#line 656 "zlahef.f"
	    i__2 = k + k * a_dim1;
#line 656 "zlahef.f"
	    d__1 = a[i__2].r;
#line 656 "zlahef.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 657 "zlahef.f"
	} else {

/*           ============================================================ */

/*           BEGIN pivot search */

/*           Case(1) */
#line 664 "zlahef.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 668 "zlahef.f"
		kp = k;
#line 669 "zlahef.f"
	    } else {

/*              BEGIN pivot search along IMAX row */


/*              Copy column IMAX to column K+1 of W and update it */

#line 676 "zlahef.f"
		i__1 = imax - k;
#line 676 "zlahef.f"
		zcopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 677 "zlahef.f"
		i__1 = imax - k;
#line 677 "zlahef.f"
		zlacgv_(&i__1, &w[k + (k + 1) * w_dim1], &c__1);
#line 678 "zlahef.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 678 "zlahef.f"
		i__2 = imax + imax * a_dim1;
#line 678 "zlahef.f"
		d__1 = a[i__2].r;
#line 678 "zlahef.f"
		w[i__1].r = d__1, w[i__1].i = 0.;
#line 679 "zlahef.f"
		if (imax < *n) {
#line 679 "zlahef.f"
		    i__1 = *n - imax;
#line 679 "zlahef.f"
		    zcopy_(&i__1, &a[imax + 1 + imax * a_dim1], &c__1, &w[
			    imax + 1 + (k + 1) * w_dim1], &c__1);
#line 679 "zlahef.f"
		}
#line 682 "zlahef.f"
		i__1 = *n - k + 1;
#line 682 "zlahef.f"
		i__2 = k - 1;
#line 682 "zlahef.f"
		z__1.r = -1., z__1.i = -0.;
#line 682 "zlahef.f"
		zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], 
			lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 1) * 
			w_dim1], &c__1, (ftnlen)12);
#line 685 "zlahef.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 685 "zlahef.f"
		i__2 = imax + (k + 1) * w_dim1;
#line 685 "zlahef.f"
		d__1 = w[i__2].r;
#line 685 "zlahef.f"
		w[i__1].r = d__1, w[i__1].i = 0.;

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value. */
/*              Determine only ROWMAX. */

#line 691 "zlahef.f"
		i__1 = imax - k;
#line 691 "zlahef.f"
		jmax = k - 1 + izamax_(&i__1, &w[k + (k + 1) * w_dim1], &c__1)
			;
#line 692 "zlahef.f"
		i__1 = jmax + (k + 1) * w_dim1;
#line 692 "zlahef.f"
		rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			jmax + (k + 1) * w_dim1]), abs(d__2));
#line 693 "zlahef.f"
		if (imax < *n) {
#line 694 "zlahef.f"
		    i__1 = *n - imax;
#line 694 "zlahef.f"
		    jmax = imax + izamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
/* Computing MAX */
#line 695 "zlahef.f"
		    i__1 = jmax + (k + 1) * w_dim1;
#line 695 "zlahef.f"
		    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&w[jmax + (k + 1) * w_dim1]), abs(
			    d__2));
#line 695 "zlahef.f"
		    rowmax = max(d__3,d__4);
#line 696 "zlahef.f"
		}

/*              Case(2) */
#line 699 "zlahef.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 703 "zlahef.f"
		    kp = k;

/*              Case(3) */
#line 706 "zlahef.f"
		} else /* if(complicated condition) */ {
#line 706 "zlahef.f"
		    i__1 = imax + (k + 1) * w_dim1;
#line 706 "zlahef.f"
		    if ((d__1 = w[i__1].r, abs(d__1)) >= alpha * rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 712 "zlahef.f"
			kp = imax;

/*                 copy column K+1 of W to column K of W */

#line 716 "zlahef.f"
			i__1 = *n - k + 1;
#line 716 "zlahef.f"
			zcopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + 
				k * w_dim1], &c__1);

/*              Case(4) */
#line 719 "zlahef.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 724 "zlahef.f"
			kp = imax;
#line 725 "zlahef.f"
			kstep = 2;
#line 726 "zlahef.f"
		    }
#line 726 "zlahef.f"
		}


/*              END pivot search along IMAX row */

#line 731 "zlahef.f"
	    }

/*           END pivot search */

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 739 "zlahef.f"
	    kk = k + kstep - 1;

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KK of W. */

#line 744 "zlahef.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K+1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 751 "zlahef.f"
		i__1 = kp + kp * a_dim1;
#line 751 "zlahef.f"
		i__2 = kk + kk * a_dim1;
#line 751 "zlahef.f"
		d__1 = a[i__2].r;
#line 751 "zlahef.f"
		a[i__1].r = d__1, a[i__1].i = 0.;
#line 752 "zlahef.f"
		i__1 = kp - kk - 1;
#line 752 "zlahef.f"
		zcopy_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 754 "zlahef.f"
		i__1 = kp - kk - 1;
#line 754 "zlahef.f"
		zlacgv_(&i__1, &a[kp + (kk + 1) * a_dim1], lda);
#line 755 "zlahef.f"
		if (kp < *n) {
#line 755 "zlahef.f"
		    i__1 = *n - kp;
#line 755 "zlahef.f"
		    zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 755 "zlahef.f"
		}

/*              Interchange rows KK and KP in first K-1 columns of A */
/*              (columns K (or K and K+1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in first KK columns of W. */

#line 763 "zlahef.f"
		if (k > 1) {
#line 763 "zlahef.f"
		    i__1 = k - 1;
#line 763 "zlahef.f"
		    zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 763 "zlahef.f"
		}
#line 765 "zlahef.f"
		zswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 766 "zlahef.f"
	    }

#line 768 "zlahef.f"
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
#line 786 "zlahef.f"
		i__1 = *n - k + 1;
#line 786 "zlahef.f"
		zcopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 787 "zlahef.f"
		if (k < *n) {

/*                 (NOTE: No need to check if A(k,k) is NOT ZERO, */
/*                  since that was ensured earlier in pivot search: */
/*                  case A(k,k) = 0 falls into 2x2 pivot case(4)) */

#line 793 "zlahef.f"
		    i__1 = k + k * a_dim1;
#line 793 "zlahef.f"
		    r1 = 1. / a[i__1].r;
#line 794 "zlahef.f"
		    i__1 = *n - k;
#line 794 "zlahef.f"
		    zdscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);

/*                 (2) Conjugate column W(k) */

#line 798 "zlahef.f"
		    i__1 = *n - k;
#line 798 "zlahef.f"
		    zlacgv_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 799 "zlahef.f"
		}

#line 801 "zlahef.f"
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

#line 818 "zlahef.f"
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

#line 862 "zlahef.f"
		    i__1 = k + 1 + k * w_dim1;
#line 862 "zlahef.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 863 "zlahef.f"
		    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
#line 863 "zlahef.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 864 "zlahef.f"
		    d_cnjg(&z__2, &d21);
#line 864 "zlahef.f"
		    z_div(&z__1, &w[k + k * w_dim1], &z__2);
#line 864 "zlahef.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 865 "zlahef.f"
		    z__1.r = d11.r * d22.r - d11.i * d22.i, z__1.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 865 "zlahef.f"
		    t = 1. / (z__1.r - 1.);
#line 866 "zlahef.f"
		    z__2.r = t, z__2.i = 0.;
#line 866 "zlahef.f"
		    z_div(&z__1, &z__2, &d21);
#line 866 "zlahef.f"
		    d21.r = z__1.r, d21.i = z__1.i;

/*                 Update elements in columns A(k) and A(k+1) as */
/*                 dot products of rows of ( W(k) W(k+1) ) and columns */
/*                 of D**(-1) */

#line 872 "zlahef.f"
		    i__1 = *n;
#line 872 "zlahef.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 873 "zlahef.f"
			i__2 = j + k * a_dim1;
#line 873 "zlahef.f"
			d_cnjg(&z__2, &d21);
#line 873 "zlahef.f"
			i__3 = j + k * w_dim1;
#line 873 "zlahef.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 873 "zlahef.f"
			i__4 = j + (k + 1) * w_dim1;
#line 873 "zlahef.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 873 "zlahef.f"
			z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
				z__2.r * z__3.i + z__2.i * z__3.r;
#line 873 "zlahef.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 875 "zlahef.f"
			i__2 = j + (k + 1) * a_dim1;
#line 875 "zlahef.f"
			i__3 = j + (k + 1) * w_dim1;
#line 875 "zlahef.f"
			z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__3.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 875 "zlahef.f"
			i__4 = j + k * w_dim1;
#line 875 "zlahef.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 875 "zlahef.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 875 "zlahef.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 876 "zlahef.f"
/* L80: */
#line 876 "zlahef.f"
		    }
#line 877 "zlahef.f"
		}

/*              Copy D(k) to A */

#line 881 "zlahef.f"
		i__1 = k + k * a_dim1;
#line 881 "zlahef.f"
		i__2 = k + k * w_dim1;
#line 881 "zlahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 882 "zlahef.f"
		i__1 = k + 1 + k * a_dim1;
#line 882 "zlahef.f"
		i__2 = k + 1 + k * w_dim1;
#line 882 "zlahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 883 "zlahef.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 883 "zlahef.f"
		i__2 = k + 1 + (k + 1) * w_dim1;
#line 883 "zlahef.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

/*              (2) Conjugate columns W(k) and W(k+1) */

#line 887 "zlahef.f"
		i__1 = *n - k;
#line 887 "zlahef.f"
		zlacgv_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 888 "zlahef.f"
		i__1 = *n - k - 1;
#line 888 "zlahef.f"
		zlacgv_(&i__1, &w[k + 2 + (k + 1) * w_dim1], &c__1);

#line 890 "zlahef.f"
	    }

#line 892 "zlahef.f"
	}

/*        Store details of the interchanges in IPIV */

#line 896 "zlahef.f"
	if (kstep == 1) {
#line 897 "zlahef.f"
	    ipiv[k] = kp;
#line 898 "zlahef.f"
	} else {
#line 899 "zlahef.f"
	    ipiv[k] = -kp;
#line 900 "zlahef.f"
	    ipiv[k + 1] = -kp;
#line 901 "zlahef.f"
	}

/*        Increase K and return to the start of the main loop */

#line 905 "zlahef.f"
	k += kstep;
#line 906 "zlahef.f"
	goto L70;

#line 908 "zlahef.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**H = A22 - L21*W**H */

/*        computing blocks of NB columns at a time (note that conjg(W) is */
/*        actually stored) */

#line 917 "zlahef.f"
	i__1 = *n;
#line 917 "zlahef.f"
	i__2 = *nb;
#line 917 "zlahef.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 918 "zlahef.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 918 "zlahef.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 922 "zlahef.f"
	    i__3 = j + jb - 1;
#line 922 "zlahef.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 923 "zlahef.f"
		i__4 = jj + jj * a_dim1;
#line 923 "zlahef.f"
		i__5 = jj + jj * a_dim1;
#line 923 "zlahef.f"
		d__1 = a[i__5].r;
#line 923 "zlahef.f"
		a[i__4].r = d__1, a[i__4].i = 0.;
#line 924 "zlahef.f"
		i__4 = j + jb - jj;
#line 924 "zlahef.f"
		i__5 = k - 1;
#line 924 "zlahef.f"
		z__1.r = -1., z__1.i = -0.;
#line 924 "zlahef.f"
		zgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1]
			, &c__1, (ftnlen)12);
#line 927 "zlahef.f"
		i__4 = jj + jj * a_dim1;
#line 927 "zlahef.f"
		i__5 = jj + jj * a_dim1;
#line 927 "zlahef.f"
		d__1 = a[i__5].r;
#line 927 "zlahef.f"
		a[i__4].r = d__1, a[i__4].i = 0.;
#line 928 "zlahef.f"
/* L100: */
#line 928 "zlahef.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 932 "zlahef.f"
	    if (j + jb <= *n) {
#line 932 "zlahef.f"
		i__3 = *n - j - jb + 1;
#line 932 "zlahef.f"
		i__4 = k - 1;
#line 932 "zlahef.f"
		z__1.r = -1., z__1.i = -0.;
#line 932 "zlahef.f"
		zgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, 
			&a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 932 "zlahef.f"
	    }
#line 936 "zlahef.f"
/* L110: */
#line 936 "zlahef.f"
	}

/*        Put L21 in standard form by partially undoing the interchanges */
/*        of rows in columns 1:k-1 looping backwards from k-1 to 1 */

#line 941 "zlahef.f"
	j = k - 1;
#line 942 "zlahef.f"
L120:

/*           Undo the interchanges (if any) of rows JJ and JP at each */
/*           step J */

/*           (Here, J is a diagonal index) */
#line 948 "zlahef.f"
	jj = j;
#line 949 "zlahef.f"
	jp = ipiv[j];
#line 950 "zlahef.f"
	if (jp < 0) {
#line 951 "zlahef.f"
	    jp = -jp;
/*              (Here, J is a diagonal index) */
#line 953 "zlahef.f"
	    --j;
#line 954 "zlahef.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length J */
/*           of the rows to swap back doesn't include diagonal element) */
#line 957 "zlahef.f"
	--j;
#line 958 "zlahef.f"
	if (jp != jj && j >= 1) {
#line 958 "zlahef.f"
	    zswap_(&j, &a[jp + a_dim1], lda, &a[jj + a_dim1], lda);
#line 958 "zlahef.f"
	}
#line 960 "zlahef.f"
	if (j > 1) {
#line 960 "zlahef.f"
	    goto L120;
#line 960 "zlahef.f"
	}

/*        Set KB to the number of columns factorized */

#line 965 "zlahef.f"
	*kb = k - 1;

#line 967 "zlahef.f"
    }
#line 968 "zlahef.f"
    return 0;

/*     End of ZLAHEF */

} /* zlahef_ */

