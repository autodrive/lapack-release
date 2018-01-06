#line 1 "clasyf.f"
/* clasyf.f -- translated by f2c (version 20100827).
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

#line 1 "clasyf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLASYF computes a partial factorization of a complex symmetric matrix using the Bunch-Kaufman d
iagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLASYF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasyf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasyf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasyf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */

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
/* > CLASYF computes a partial factorization of a complex symmetric matrix */
/* > A using the Bunch-Kaufman diagonal pivoting method. The partial */
/* > factorization has the form: */
/* > */
/* > A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or: */
/* >       ( 0  U22 ) (  0   D  ) ( U12**T U22**T ) */
/* > */
/* > A  =  ( L11  0 ) ( D    0  ) ( L11**T L21**T )  if UPLO = 'L' */
/* >       ( L21  I ) ( 0   A22 ) (  0       I    ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > Note that U**T denotes the transpose of U. */
/* > */
/* > CLASYF is an auxiliary routine called by CSYTRF. It uses blocked code */
/* > (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or */
/* > A22 (if UPLO = 'L'). */
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
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
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

/* > \ingroup complexSYcomputational */

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
/* Subroutine */ int clasyf_(char *uplo, integer *n, integer *nb, integer *kb,
	 doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, 
	integer *ldw, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex t, r1, d11, d21, d22;
    static integer jb, jj, kk, jp, kp, kw, kkw, imax, jmax;
    static doublereal alpha;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
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
    extern integer icamax_(integer *, doublecomplex *, integer *);
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

#line 229 "clasyf.f"
    /* Parameter adjustments */
#line 229 "clasyf.f"
    a_dim1 = *lda;
#line 229 "clasyf.f"
    a_offset = 1 + a_dim1;
#line 229 "clasyf.f"
    a -= a_offset;
#line 229 "clasyf.f"
    --ipiv;
#line 229 "clasyf.f"
    w_dim1 = *ldw;
#line 229 "clasyf.f"
    w_offset = 1 + w_dim1;
#line 229 "clasyf.f"
    w -= w_offset;
#line 229 "clasyf.f"

#line 229 "clasyf.f"
    /* Function Body */
#line 229 "clasyf.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 233 "clasyf.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 235 "clasyf.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 */

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

/*        KW is the column of W which corresponds to column K of A */

#line 245 "clasyf.f"
	k = *n;
#line 246 "clasyf.f"
L10:
#line 247 "clasyf.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 251 "clasyf.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 251 "clasyf.f"
	    goto L30;
#line 251 "clasyf.f"
	}

/*        Copy column K of A to column KW of W and update it */

#line 256 "clasyf.f"
	ccopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
#line 257 "clasyf.f"
	if (k < *n) {
#line 257 "clasyf.f"
	    i__1 = *n - k;
#line 257 "clasyf.f"
	    z__1.r = -1., z__1.i = -0.;
#line 257 "clasyf.f"
	    cgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 257 "clasyf.f"
	}

#line 261 "clasyf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 266 "clasyf.f"
	i__1 = k + kw * w_dim1;
#line 266 "clasyf.f"
	absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[k + kw * 
		w_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 272 "clasyf.f"
	if (k > 1) {
#line 273 "clasyf.f"
	    i__1 = k - 1;
#line 273 "clasyf.f"
	    imax = icamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 274 "clasyf.f"
	    i__1 = imax + kw * w_dim1;
#line 274 "clasyf.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    kw * w_dim1]), abs(d__2));
#line 275 "clasyf.f"
	} else {
#line 276 "clasyf.f"
	    colmax = 0.;
#line 277 "clasyf.f"
	}

#line 279 "clasyf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 283 "clasyf.f"
	    if (*info == 0) {
#line 283 "clasyf.f"
		*info = k;
#line 283 "clasyf.f"
	    }
#line 285 "clasyf.f"
	    kp = k;
#line 286 "clasyf.f"
	} else {
#line 287 "clasyf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 291 "clasyf.f"
		kp = k;
#line 292 "clasyf.f"
	    } else {

/*              Copy column IMAX to column KW-1 of W and update it */

#line 296 "clasyf.f"
		ccopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
#line 297 "clasyf.f"
		i__1 = k - imax;
#line 297 "clasyf.f"
		ccopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);
#line 299 "clasyf.f"
		if (k < *n) {
#line 299 "clasyf.f"
		    i__1 = *n - k;
#line 299 "clasyf.f"
		    z__1.r = -1., z__1.i = -0.;
#line 299 "clasyf.f"
		    cgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 299 "clasyf.f"
		}

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 307 "clasyf.f"
		i__1 = k - imax;
#line 307 "clasyf.f"
		jmax = imax + icamax_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1],
			 &c__1);
#line 308 "clasyf.f"
		i__1 = jmax + (kw - 1) * w_dim1;
#line 308 "clasyf.f"
		rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			jmax + (kw - 1) * w_dim1]), abs(d__2));
#line 309 "clasyf.f"
		if (imax > 1) {
#line 310 "clasyf.f"
		    i__1 = imax - 1;
#line 310 "clasyf.f"
		    jmax = icamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
/* Computing MAX */
#line 311 "clasyf.f"
		    i__1 = jmax + (kw - 1) * w_dim1;
#line 311 "clasyf.f"
		    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&w[jmax + (kw - 1) * w_dim1]), abs(
			    d__2));
#line 311 "clasyf.f"
		    rowmax = max(d__3,d__4);
#line 312 "clasyf.f"
		}

#line 314 "clasyf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 318 "clasyf.f"
		    kp = k;
#line 319 "clasyf.f"
		} else /* if(complicated condition) */ {
#line 319 "clasyf.f"
		    i__1 = imax + (kw - 1) * w_dim1;
#line 319 "clasyf.f"
		    if ((d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    imax + (kw - 1) * w_dim1]), abs(d__2)) >= alpha * 
			    rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 324 "clasyf.f"
			kp = imax;

/*                 copy column KW-1 of W to column KW of W */

#line 328 "clasyf.f"
			ccopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
				w_dim1 + 1], &c__1);
#line 329 "clasyf.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 334 "clasyf.f"
			kp = imax;
#line 335 "clasyf.f"
			kstep = 2;
#line 336 "clasyf.f"
		    }
#line 336 "clasyf.f"
		}
#line 337 "clasyf.f"
	    }

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 343 "clasyf.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 347 "clasyf.f"
	    kkw = *nb + kk - *n;

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KKW of W. */

#line 352 "clasyf.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K-1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 359 "clasyf.f"
		i__1 = kp + kp * a_dim1;
#line 359 "clasyf.f"
		i__2 = kk + kk * a_dim1;
#line 359 "clasyf.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 360 "clasyf.f"
		i__1 = kk - 1 - kp;
#line 360 "clasyf.f"
		ccopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 362 "clasyf.f"
		if (kp > 1) {
#line 362 "clasyf.f"
		    i__1 = kp - 1;
#line 362 "clasyf.f"
		    ccopy_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 362 "clasyf.f"
		}

/*              Interchange rows KK and KP in last K+1 to N columns of A */
/*              (columns K (or K and K-1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in last KKW to NB columns of W. */

#line 370 "clasyf.f"
		if (k < *n) {
#line 370 "clasyf.f"
		    i__1 = *n - k;
#line 370 "clasyf.f"
		    cswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 370 "clasyf.f"
		}
#line 373 "clasyf.f"
		i__1 = *n - kk + 1;
#line 373 "clasyf.f"
		cswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 375 "clasyf.f"
	    }

#line 377 "clasyf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column kw of W now holds */

/*              W(kw) = U(k)*D(k), */

/*              where U(k) is the k-th column of U */

/*              Store subdiag. elements of column U(k) */
/*              and 1-by-1 block D(k) in column k of A. */
/*              NOTE: Diagonal element U(k,k) is a UNIT element */
/*              and not stored. */
/*                 A(k,k) := D(k,k) = W(k,kw) */
/*                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k) */

#line 392 "clasyf.f"
		ccopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 393 "clasyf.f"
		z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 393 "clasyf.f"
		r1.r = z__1.r, r1.i = z__1.i;
#line 394 "clasyf.f"
		i__1 = k - 1;
#line 394 "clasyf.f"
		cscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);

#line 396 "clasyf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns kw and kw-1 of W now hold */

/*              ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

/*              Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2 */
/*              block D(k-1:k,k-1:k) in columns k-1 and k of A. */
/*              NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT */
/*              block and not stored. */
/*                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw) */
/*                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) = */
/*                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) ) */

#line 413 "clasyf.f"
		if (k > 2) {

/*                 Compose the columns of the inverse of 2-by-2 pivot */
/*                 block D in the following way to reduce the number */
/*                 of FLOPS when we myltiply panel ( W(kw-1) W(kw) ) by */
/*                 this inverse */

/*                 D**(-1) = ( d11 d21 )**(-1) = */
/*                           ( d21 d22 ) */

/*                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) = */
/*                                        ( (-d21 ) ( d11 ) ) */

/*                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) * */

/*                   * ( ( d22/d21 ) (      -1 ) ) = */
/*                     ( (      -1 ) ( d11/d21 ) ) */

/*                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) = */
/*                                           ( ( -1  ) ( D22 ) ) */

/*                 = 1/d21 * T * ( ( D11 ) (  -1 ) ) */
/*                               ( (  -1 ) ( D22 ) ) */

/*                 = D21 * ( ( D11 ) (  -1 ) ) */
/*                         ( (  -1 ) ( D22 ) ) */

#line 440 "clasyf.f"
		    i__1 = k - 1 + kw * w_dim1;
#line 440 "clasyf.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 441 "clasyf.f"
		    z_div(&z__1, &w[k + kw * w_dim1], &d21);
#line 441 "clasyf.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 442 "clasyf.f"
		    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d21);
#line 442 "clasyf.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 443 "clasyf.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 443 "clasyf.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 443 "clasyf.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 443 "clasyf.f"
		    t.r = z__1.r, t.i = z__1.i;

/*                 Update elements in columns A(k-1) and A(k) as */
/*                 dot products of rows of ( W(kw-1) W(kw) ) and columns */
/*                 of D**(-1) */

#line 449 "clasyf.f"
		    z_div(&z__1, &t, &d21);
#line 449 "clasyf.f"
		    d21.r = z__1.r, d21.i = z__1.i;
#line 450 "clasyf.f"
		    i__1 = k - 2;
#line 450 "clasyf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 451 "clasyf.f"
			i__2 = j + (k - 1) * a_dim1;
#line 451 "clasyf.f"
			i__3 = j + (kw - 1) * w_dim1;
#line 451 "clasyf.f"
			z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__3.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 451 "clasyf.f"
			i__4 = j + kw * w_dim1;
#line 451 "clasyf.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 451 "clasyf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 451 "clasyf.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 452 "clasyf.f"
			i__2 = j + k * a_dim1;
#line 452 "clasyf.f"
			i__3 = j + kw * w_dim1;
#line 452 "clasyf.f"
			z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__3.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 452 "clasyf.f"
			i__4 = j + (kw - 1) * w_dim1;
#line 452 "clasyf.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 452 "clasyf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 452 "clasyf.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 453 "clasyf.f"
/* L20: */
#line 453 "clasyf.f"
		    }
#line 454 "clasyf.f"
		}

/*              Copy D(k) to A */

#line 458 "clasyf.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 458 "clasyf.f"
		i__2 = k - 1 + (kw - 1) * w_dim1;
#line 458 "clasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 459 "clasyf.f"
		i__1 = k - 1 + k * a_dim1;
#line 459 "clasyf.f"
		i__2 = k - 1 + kw * w_dim1;
#line 459 "clasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 460 "clasyf.f"
		i__1 = k + k * a_dim1;
#line 460 "clasyf.f"
		i__2 = k + kw * w_dim1;
#line 460 "clasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

#line 462 "clasyf.f"
	    }

#line 464 "clasyf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 468 "clasyf.f"
	if (kstep == 1) {
#line 469 "clasyf.f"
	    ipiv[k] = kp;
#line 470 "clasyf.f"
	} else {
#line 471 "clasyf.f"
	    ipiv[k] = -kp;
#line 472 "clasyf.f"
	    ipiv[k - 1] = -kp;
#line 473 "clasyf.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 477 "clasyf.f"
	k -= kstep;
#line 478 "clasyf.f"
	goto L10;

#line 480 "clasyf.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T */

/*        computing blocks of NB columns at a time */

#line 488 "clasyf.f"
	i__1 = -(*nb);
#line 488 "clasyf.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 489 "clasyf.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 489 "clasyf.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 493 "clasyf.f"
	    i__2 = j + jb - 1;
#line 493 "clasyf.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 494 "clasyf.f"
		i__3 = jj - j + 1;
#line 494 "clasyf.f"
		i__4 = *n - k;
#line 494 "clasyf.f"
		z__1.r = -1., z__1.i = -0.;
#line 494 "clasyf.f"
		cgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, 
			&a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 497 "clasyf.f"
/* L40: */
#line 497 "clasyf.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 501 "clasyf.f"
	    i__2 = j - 1;
#line 501 "clasyf.f"
	    i__3 = *n - k;
#line 501 "clasyf.f"
	    z__1.r = -1., z__1.i = -0.;
#line 501 "clasyf.f"
	    cgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, &a[(
		    k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * w_dim1], ldw,
		     &c_b1, &a[j * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
#line 504 "clasyf.f"
/* L50: */
#line 504 "clasyf.f"
	}

/*        Put U12 in standard form by partially undoing the interchanges */
/*        in columns k+1:n looping backwards from k+1 to n */

#line 509 "clasyf.f"
	j = k + 1;
#line 510 "clasyf.f"
L60:

/*           Undo the interchanges (if any) of rows JJ and JP at each */
/*           step J */

/*           (Here, J is a diagonal index) */
#line 516 "clasyf.f"
	jj = j;
#line 517 "clasyf.f"
	jp = ipiv[j];
#line 518 "clasyf.f"
	if (jp < 0) {
#line 519 "clasyf.f"
	    jp = -jp;
/*              (Here, J is a diagonal index) */
#line 521 "clasyf.f"
	    ++j;
#line 522 "clasyf.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length N-J+1 */
/*           of the rows to swap back doesn't include diagonal element) */
#line 525 "clasyf.f"
	++j;
#line 526 "clasyf.f"
	if (jp != jj && j <= *n) {
#line 526 "clasyf.f"
	    i__1 = *n - j + 1;
#line 526 "clasyf.f"
	    cswap_(&i__1, &a[jp + j * a_dim1], lda, &a[jj + j * a_dim1], lda);
#line 526 "clasyf.f"
	}
#line 528 "clasyf.f"
	if (j < *n) {
#line 528 "clasyf.f"
	    goto L60;
#line 528 "clasyf.f"
	}

/*        Set KB to the number of columns factorized */

#line 533 "clasyf.f"
	*kb = *n - k;

#line 535 "clasyf.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 */

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 543 "clasyf.f"
	k = 1;
#line 544 "clasyf.f"
L70:

/*        Exit from loop */

#line 548 "clasyf.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 548 "clasyf.f"
	    goto L90;
#line 548 "clasyf.f"
	}

/*        Copy column K of A to column K of W and update it */

#line 553 "clasyf.f"
	i__1 = *n - k + 1;
#line 553 "clasyf.f"
	ccopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
#line 554 "clasyf.f"
	i__1 = *n - k + 1;
#line 554 "clasyf.f"
	i__2 = k - 1;
#line 554 "clasyf.f"
	z__1.r = -1., z__1.i = -0.;
#line 554 "clasyf.f"
	cgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[k 
		+ w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1, (ftnlen)12);

#line 557 "clasyf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 562 "clasyf.f"
	i__1 = k + k * w_dim1;
#line 562 "clasyf.f"
	absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[k + k * 
		w_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 568 "clasyf.f"
	if (k < *n) {
#line 569 "clasyf.f"
	    i__1 = *n - k;
#line 569 "clasyf.f"
	    imax = k + icamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 570 "clasyf.f"
	    i__1 = imax + k * w_dim1;
#line 570 "clasyf.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    k * w_dim1]), abs(d__2));
#line 571 "clasyf.f"
	} else {
#line 572 "clasyf.f"
	    colmax = 0.;
#line 573 "clasyf.f"
	}

#line 575 "clasyf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 579 "clasyf.f"
	    if (*info == 0) {
#line 579 "clasyf.f"
		*info = k;
#line 579 "clasyf.f"
	    }
#line 581 "clasyf.f"
	    kp = k;
#line 582 "clasyf.f"
	} else {
#line 583 "clasyf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 587 "clasyf.f"
		kp = k;
#line 588 "clasyf.f"
	    } else {

/*              Copy column IMAX to column K+1 of W and update it */

#line 592 "clasyf.f"
		i__1 = imax - k;
#line 592 "clasyf.f"
		ccopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 593 "clasyf.f"
		i__1 = *n - imax + 1;
#line 593 "clasyf.f"
		ccopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 
			1) * w_dim1], &c__1);
#line 595 "clasyf.f"
		i__1 = *n - k + 1;
#line 595 "clasyf.f"
		i__2 = k - 1;
#line 595 "clasyf.f"
		z__1.r = -1., z__1.i = -0.;
#line 595 "clasyf.f"
		cgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], 
			lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 1) * 
			w_dim1], &c__1, (ftnlen)12);

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 602 "clasyf.f"
		i__1 = imax - k;
#line 602 "clasyf.f"
		jmax = k - 1 + icamax_(&i__1, &w[k + (k + 1) * w_dim1], &c__1)
			;
#line 603 "clasyf.f"
		i__1 = jmax + (k + 1) * w_dim1;
#line 603 "clasyf.f"
		rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			jmax + (k + 1) * w_dim1]), abs(d__2));
#line 604 "clasyf.f"
		if (imax < *n) {
#line 605 "clasyf.f"
		    i__1 = *n - imax;
#line 605 "clasyf.f"
		    jmax = imax + icamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
/* Computing MAX */
#line 606 "clasyf.f"
		    i__1 = jmax + (k + 1) * w_dim1;
#line 606 "clasyf.f"
		    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&w[jmax + (k + 1) * w_dim1]), abs(
			    d__2));
#line 606 "clasyf.f"
		    rowmax = max(d__3,d__4);
#line 607 "clasyf.f"
		}

#line 609 "clasyf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 613 "clasyf.f"
		    kp = k;
#line 614 "clasyf.f"
		} else /* if(complicated condition) */ {
#line 614 "clasyf.f"
		    i__1 = imax + (k + 1) * w_dim1;
#line 614 "clasyf.f"
		    if ((d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    imax + (k + 1) * w_dim1]), abs(d__2)) >= alpha * 
			    rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 619 "clasyf.f"
			kp = imax;

/*                 copy column K+1 of W to column K of W */

#line 623 "clasyf.f"
			i__1 = *n - k + 1;
#line 623 "clasyf.f"
			ccopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + 
				k * w_dim1], &c__1);
#line 624 "clasyf.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 629 "clasyf.f"
			kp = imax;
#line 630 "clasyf.f"
			kstep = 2;
#line 631 "clasyf.f"
		    }
#line 631 "clasyf.f"
		}
#line 632 "clasyf.f"
	    }

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 638 "clasyf.f"
	    kk = k + kstep - 1;

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KK of W. */

#line 643 "clasyf.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K+1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 650 "clasyf.f"
		i__1 = kp + kp * a_dim1;
#line 650 "clasyf.f"
		i__2 = kk + kk * a_dim1;
#line 650 "clasyf.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 651 "clasyf.f"
		i__1 = kp - kk - 1;
#line 651 "clasyf.f"
		ccopy_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 653 "clasyf.f"
		if (kp < *n) {
#line 653 "clasyf.f"
		    i__1 = *n - kp;
#line 653 "clasyf.f"
		    ccopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 653 "clasyf.f"
		}

/*              Interchange rows KK and KP in first K-1 columns of A */
/*              (columns K (or K and K+1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in first KK columns of W. */

#line 661 "clasyf.f"
		if (k > 1) {
#line 661 "clasyf.f"
		    i__1 = k - 1;
#line 661 "clasyf.f"
		    cswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 661 "clasyf.f"
		}
#line 663 "clasyf.f"
		cswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 664 "clasyf.f"
	    }

#line 666 "clasyf.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of W now holds */

/*              W(k) = L(k)*D(k), */

/*              where L(k) is the k-th column of L */

/*              Store subdiag. elements of column L(k) */
/*              and 1-by-1 block D(k) in column k of A. */
/*              (NOTE: Diagonal element L(k,k) is a UNIT element */
/*              and not stored) */
/*                 A(k,k) := D(k,k) = W(k,k) */
/*                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k) */

#line 681 "clasyf.f"
		i__1 = *n - k + 1;
#line 681 "clasyf.f"
		ccopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 682 "clasyf.f"
		if (k < *n) {
#line 683 "clasyf.f"
		    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 683 "clasyf.f"
		    r1.r = z__1.r, r1.i = z__1.i;
#line 684 "clasyf.f"
		    i__1 = *n - k;
#line 684 "clasyf.f"
		    cscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 685 "clasyf.f"
		}

#line 687 "clasyf.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

/*              Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2 */
/*              block D(k:k+1,k:k+1) in columns k and k+1 of A. */
/*              (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT */
/*              block and not stored) */
/*                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1) */
/*                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) = */
/*                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) ) */

#line 704 "clasyf.f"
		if (k < *n - 1) {

/*                 Compose the columns of the inverse of 2-by-2 pivot */
/*                 block D in the following way to reduce the number */
/*                 of FLOPS when we myltiply panel ( W(k) W(k+1) ) by */
/*                 this inverse */

/*                 D**(-1) = ( d11 d21 )**(-1) = */
/*                           ( d21 d22 ) */

/*                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) = */
/*                                        ( (-d21 ) ( d11 ) ) */

/*                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) * */

/*                   * ( ( d22/d21 ) (      -1 ) ) = */
/*                     ( (      -1 ) ( d11/d21 ) ) */

/*                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) = */
/*                                           ( ( -1  ) ( D22 ) ) */

/*                 = 1/d21 * T * ( ( D11 ) (  -1 ) ) */
/*                               ( (  -1 ) ( D22 ) ) */

/*                 = D21 * ( ( D11 ) (  -1 ) ) */
/*                         ( (  -1 ) ( D22 ) ) */

#line 731 "clasyf.f"
		    i__1 = k + 1 + k * w_dim1;
#line 731 "clasyf.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 732 "clasyf.f"
		    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
#line 732 "clasyf.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 733 "clasyf.f"
		    z_div(&z__1, &w[k + k * w_dim1], &d21);
#line 733 "clasyf.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 734 "clasyf.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 734 "clasyf.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 734 "clasyf.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 734 "clasyf.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 735 "clasyf.f"
		    z_div(&z__1, &t, &d21);
#line 735 "clasyf.f"
		    d21.r = z__1.r, d21.i = z__1.i;

/*                 Update elements in columns A(k) and A(k+1) as */
/*                 dot products of rows of ( W(k) W(k+1) ) and columns */
/*                 of D**(-1) */

#line 741 "clasyf.f"
		    i__1 = *n;
#line 741 "clasyf.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 742 "clasyf.f"
			i__2 = j + k * a_dim1;
#line 742 "clasyf.f"
			i__3 = j + k * w_dim1;
#line 742 "clasyf.f"
			z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__3.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 742 "clasyf.f"
			i__4 = j + (k + 1) * w_dim1;
#line 742 "clasyf.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 742 "clasyf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 742 "clasyf.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 743 "clasyf.f"
			i__2 = j + (k + 1) * a_dim1;
#line 743 "clasyf.f"
			i__3 = j + (k + 1) * w_dim1;
#line 743 "clasyf.f"
			z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__3.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 743 "clasyf.f"
			i__4 = j + k * w_dim1;
#line 743 "clasyf.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 743 "clasyf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 743 "clasyf.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 744 "clasyf.f"
/* L80: */
#line 744 "clasyf.f"
		    }
#line 745 "clasyf.f"
		}

/*              Copy D(k) to A */

#line 749 "clasyf.f"
		i__1 = k + k * a_dim1;
#line 749 "clasyf.f"
		i__2 = k + k * w_dim1;
#line 749 "clasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 750 "clasyf.f"
		i__1 = k + 1 + k * a_dim1;
#line 750 "clasyf.f"
		i__2 = k + 1 + k * w_dim1;
#line 750 "clasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 751 "clasyf.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 751 "clasyf.f"
		i__2 = k + 1 + (k + 1) * w_dim1;
#line 751 "clasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

#line 753 "clasyf.f"
	    }

#line 755 "clasyf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 759 "clasyf.f"
	if (kstep == 1) {
#line 760 "clasyf.f"
	    ipiv[k] = kp;
#line 761 "clasyf.f"
	} else {
#line 762 "clasyf.f"
	    ipiv[k] = -kp;
#line 763 "clasyf.f"
	    ipiv[k + 1] = -kp;
#line 764 "clasyf.f"
	}

/*        Increase K and return to the start of the main loop */

#line 768 "clasyf.f"
	k += kstep;
#line 769 "clasyf.f"
	goto L70;

#line 771 "clasyf.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T */

/*        computing blocks of NB columns at a time */

#line 779 "clasyf.f"
	i__1 = *n;
#line 779 "clasyf.f"
	i__2 = *nb;
#line 779 "clasyf.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 780 "clasyf.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 780 "clasyf.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 784 "clasyf.f"
	    i__3 = j + jb - 1;
#line 784 "clasyf.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 785 "clasyf.f"
		i__4 = j + jb - jj;
#line 785 "clasyf.f"
		i__5 = k - 1;
#line 785 "clasyf.f"
		z__1.r = -1., z__1.i = -0.;
#line 785 "clasyf.f"
		cgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1]
			, &c__1, (ftnlen)12);
#line 788 "clasyf.f"
/* L100: */
#line 788 "clasyf.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 792 "clasyf.f"
	    if (j + jb <= *n) {
#line 792 "clasyf.f"
		i__3 = *n - j - jb + 1;
#line 792 "clasyf.f"
		i__4 = k - 1;
#line 792 "clasyf.f"
		z__1.r = -1., z__1.i = -0.;
#line 792 "clasyf.f"
		cgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, 
			&a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 792 "clasyf.f"
	    }
#line 796 "clasyf.f"
/* L110: */
#line 796 "clasyf.f"
	}

/*        Put L21 in standard form by partially undoing the interchanges */
/*        of rows in columns 1:k-1 looping backwards from k-1 to 1 */

#line 801 "clasyf.f"
	j = k - 1;
#line 802 "clasyf.f"
L120:

/*           Undo the interchanges (if any) of rows JJ and JP at each */
/*           step J */

/*           (Here, J is a diagonal index) */
#line 808 "clasyf.f"
	jj = j;
#line 809 "clasyf.f"
	jp = ipiv[j];
#line 810 "clasyf.f"
	if (jp < 0) {
#line 811 "clasyf.f"
	    jp = -jp;
/*              (Here, J is a diagonal index) */
#line 813 "clasyf.f"
	    --j;
#line 814 "clasyf.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length J */
/*           of the rows to swap back doesn't include diagonal element) */
#line 817 "clasyf.f"
	--j;
#line 818 "clasyf.f"
	if (jp != jj && j >= 1) {
#line 818 "clasyf.f"
	    cswap_(&j, &a[jp + a_dim1], lda, &a[jj + a_dim1], lda);
#line 818 "clasyf.f"
	}
#line 820 "clasyf.f"
	if (j > 1) {
#line 820 "clasyf.f"
	    goto L120;
#line 820 "clasyf.f"
	}

/*        Set KB to the number of columns factorized */

#line 825 "clasyf.f"
	*kb = k - 1;

#line 827 "clasyf.f"
    }
#line 828 "clasyf.f"
    return 0;

/*     End of CLASYF */

} /* clasyf_ */

