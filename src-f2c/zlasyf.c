#line 1 "zlasyf.f"
/* zlasyf.f -- translated by f2c (version 20100827).
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

#line 1 "zlasyf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLASYF computes a partial factorization of a complex symmetric matrix using the Bunch-Kaufman d
iagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLASYF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlasyf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlasyf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlasyf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */

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
/* > ZLASYF computes a partial factorization of a complex symmetric matrix */
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
/* > ZLASYF is an auxiliary routine called by ZSYTRF. It uses blocked code */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \ingroup complex16SYcomputational */

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
/* Subroutine */ int zlasyf_(char *uplo, integer *n, integer *nb, integer *kb,
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer kstep;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal absakk, colmax;
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

#line 229 "zlasyf.f"
    /* Parameter adjustments */
#line 229 "zlasyf.f"
    a_dim1 = *lda;
#line 229 "zlasyf.f"
    a_offset = 1 + a_dim1;
#line 229 "zlasyf.f"
    a -= a_offset;
#line 229 "zlasyf.f"
    --ipiv;
#line 229 "zlasyf.f"
    w_dim1 = *ldw;
#line 229 "zlasyf.f"
    w_offset = 1 + w_dim1;
#line 229 "zlasyf.f"
    w -= w_offset;
#line 229 "zlasyf.f"

#line 229 "zlasyf.f"
    /* Function Body */
#line 229 "zlasyf.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 233 "zlasyf.f"
    alpha = (sqrt(17.) + 1.) / 8.;

#line 235 "zlasyf.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 */

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

/*        KW is the column of W which corresponds to column K of A */

#line 245 "zlasyf.f"
	k = *n;
#line 246 "zlasyf.f"
L10:
#line 247 "zlasyf.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 251 "zlasyf.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 251 "zlasyf.f"
	    goto L30;
#line 251 "zlasyf.f"
	}

/*        Copy column K of A to column KW of W and update it */

#line 256 "zlasyf.f"
	zcopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
#line 257 "zlasyf.f"
	if (k < *n) {
#line 257 "zlasyf.f"
	    i__1 = *n - k;
#line 257 "zlasyf.f"
	    z__1.r = -1., z__1.i = -0.;
#line 257 "zlasyf.f"
	    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 257 "zlasyf.f"
	}

#line 261 "zlasyf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 266 "zlasyf.f"
	i__1 = k + kw * w_dim1;
#line 266 "zlasyf.f"
	absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[k + kw * 
		w_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */

#line 271 "zlasyf.f"
	if (k > 1) {
#line 272 "zlasyf.f"
	    i__1 = k - 1;
#line 272 "zlasyf.f"
	    imax = izamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 273 "zlasyf.f"
	    i__1 = imax + kw * w_dim1;
#line 273 "zlasyf.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    kw * w_dim1]), abs(d__2));
#line 274 "zlasyf.f"
	} else {
#line 275 "zlasyf.f"
	    colmax = 0.;
#line 276 "zlasyf.f"
	}

#line 278 "zlasyf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 282 "zlasyf.f"
	    if (*info == 0) {
#line 282 "zlasyf.f"
		*info = k;
#line 282 "zlasyf.f"
	    }
#line 284 "zlasyf.f"
	    kp = k;
#line 285 "zlasyf.f"
	} else {
#line 286 "zlasyf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 290 "zlasyf.f"
		kp = k;
#line 291 "zlasyf.f"
	    } else {

/*              Copy column IMAX to column KW-1 of W and update it */

#line 295 "zlasyf.f"
		zcopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
#line 296 "zlasyf.f"
		i__1 = k - imax;
#line 296 "zlasyf.f"
		zcopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);
#line 298 "zlasyf.f"
		if (k < *n) {
#line 298 "zlasyf.f"
		    i__1 = *n - k;
#line 298 "zlasyf.f"
		    z__1.r = -1., z__1.i = -0.;
#line 298 "zlasyf.f"
		    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 298 "zlasyf.f"
		}

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 306 "zlasyf.f"
		i__1 = k - imax;
#line 306 "zlasyf.f"
		jmax = imax + izamax_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1],
			 &c__1);
#line 307 "zlasyf.f"
		i__1 = jmax + (kw - 1) * w_dim1;
#line 307 "zlasyf.f"
		rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			jmax + (kw - 1) * w_dim1]), abs(d__2));
#line 308 "zlasyf.f"
		if (imax > 1) {
#line 309 "zlasyf.f"
		    i__1 = imax - 1;
#line 309 "zlasyf.f"
		    jmax = izamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
/* Computing MAX */
#line 310 "zlasyf.f"
		    i__1 = jmax + (kw - 1) * w_dim1;
#line 310 "zlasyf.f"
		    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&w[jmax + (kw - 1) * w_dim1]), abs(
			    d__2));
#line 310 "zlasyf.f"
		    rowmax = max(d__3,d__4);
#line 311 "zlasyf.f"
		}

#line 313 "zlasyf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 317 "zlasyf.f"
		    kp = k;
#line 318 "zlasyf.f"
		} else /* if(complicated condition) */ {
#line 318 "zlasyf.f"
		    i__1 = imax + (kw - 1) * w_dim1;
#line 318 "zlasyf.f"
		    if ((d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    imax + (kw - 1) * w_dim1]), abs(d__2)) >= alpha * 
			    rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 323 "zlasyf.f"
			kp = imax;

/*                 copy column KW-1 of W to column KW of W */

#line 327 "zlasyf.f"
			zcopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
				w_dim1 + 1], &c__1);
#line 328 "zlasyf.f"
		    } else {

/*                 interchange rows and columns K-1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 333 "zlasyf.f"
			kp = imax;
#line 334 "zlasyf.f"
			kstep = 2;
#line 335 "zlasyf.f"
		    }
#line 335 "zlasyf.f"
		}
#line 336 "zlasyf.f"
	    }

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 342 "zlasyf.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 346 "zlasyf.f"
	    kkw = *nb + kk - *n;

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KKW of W. */

#line 351 "zlasyf.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K-1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 358 "zlasyf.f"
		i__1 = kp + kp * a_dim1;
#line 358 "zlasyf.f"
		i__2 = kk + kk * a_dim1;
#line 358 "zlasyf.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 359 "zlasyf.f"
		i__1 = kk - 1 - kp;
#line 359 "zlasyf.f"
		zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 361 "zlasyf.f"
		if (kp > 1) {
#line 361 "zlasyf.f"
		    i__1 = kp - 1;
#line 361 "zlasyf.f"
		    zcopy_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 
			    + 1], &c__1);
#line 361 "zlasyf.f"
		}

/*              Interchange rows KK and KP in last K+1 to N columns of A */
/*              (columns K (or K and K-1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in last KKW to NB columns of W. */

#line 369 "zlasyf.f"
		if (k < *n) {
#line 369 "zlasyf.f"
		    i__1 = *n - k;
#line 369 "zlasyf.f"
		    zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k 
			    + 1) * a_dim1], lda);
#line 369 "zlasyf.f"
		}
#line 372 "zlasyf.f"
		i__1 = *n - kk + 1;
#line 372 "zlasyf.f"
		zswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 374 "zlasyf.f"
	    }

#line 376 "zlasyf.f"
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

#line 391 "zlasyf.f"
		zcopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 392 "zlasyf.f"
		z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 392 "zlasyf.f"
		r1.r = z__1.r, r1.i = z__1.i;
#line 393 "zlasyf.f"
		i__1 = k - 1;
#line 393 "zlasyf.f"
		zscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);

#line 395 "zlasyf.f"
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

#line 412 "zlasyf.f"
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

#line 439 "zlasyf.f"
		    i__1 = k - 1 + kw * w_dim1;
#line 439 "zlasyf.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 440 "zlasyf.f"
		    z_div(&z__1, &w[k + kw * w_dim1], &d21);
#line 440 "zlasyf.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 441 "zlasyf.f"
		    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d21);
#line 441 "zlasyf.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 442 "zlasyf.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 442 "zlasyf.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 442 "zlasyf.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 442 "zlasyf.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 443 "zlasyf.f"
		    z_div(&z__1, &t, &d21);
#line 443 "zlasyf.f"
		    d21.r = z__1.r, d21.i = z__1.i;

/*                 Update elements in columns A(k-1) and A(k) as */
/*                 dot products of rows of ( W(kw-1) W(kw) ) and columns */
/*                 of D**(-1) */

#line 449 "zlasyf.f"
		    i__1 = k - 2;
#line 449 "zlasyf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 450 "zlasyf.f"
			i__2 = j + (k - 1) * a_dim1;
#line 450 "zlasyf.f"
			i__3 = j + (kw - 1) * w_dim1;
#line 450 "zlasyf.f"
			z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__3.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 450 "zlasyf.f"
			i__4 = j + kw * w_dim1;
#line 450 "zlasyf.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 450 "zlasyf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 450 "zlasyf.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 451 "zlasyf.f"
			i__2 = j + k * a_dim1;
#line 451 "zlasyf.f"
			i__3 = j + kw * w_dim1;
#line 451 "zlasyf.f"
			z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__3.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 451 "zlasyf.f"
			i__4 = j + (kw - 1) * w_dim1;
#line 451 "zlasyf.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 451 "zlasyf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 451 "zlasyf.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 452 "zlasyf.f"
/* L20: */
#line 452 "zlasyf.f"
		    }
#line 453 "zlasyf.f"
		}

/*              Copy D(k) to A */

#line 457 "zlasyf.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 457 "zlasyf.f"
		i__2 = k - 1 + (kw - 1) * w_dim1;
#line 457 "zlasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 458 "zlasyf.f"
		i__1 = k - 1 + k * a_dim1;
#line 458 "zlasyf.f"
		i__2 = k - 1 + kw * w_dim1;
#line 458 "zlasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 459 "zlasyf.f"
		i__1 = k + k * a_dim1;
#line 459 "zlasyf.f"
		i__2 = k + kw * w_dim1;
#line 459 "zlasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

#line 461 "zlasyf.f"
	    }

#line 463 "zlasyf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 467 "zlasyf.f"
	if (kstep == 1) {
#line 468 "zlasyf.f"
	    ipiv[k] = kp;
#line 469 "zlasyf.f"
	} else {
#line 470 "zlasyf.f"
	    ipiv[k] = -kp;
#line 471 "zlasyf.f"
	    ipiv[k - 1] = -kp;
#line 472 "zlasyf.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 476 "zlasyf.f"
	k -= kstep;
#line 477 "zlasyf.f"
	goto L10;

#line 479 "zlasyf.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T */

/*        computing blocks of NB columns at a time */

#line 487 "zlasyf.f"
	i__1 = -(*nb);
#line 487 "zlasyf.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 488 "zlasyf.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 488 "zlasyf.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 492 "zlasyf.f"
	    i__2 = j + jb - 1;
#line 492 "zlasyf.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 493 "zlasyf.f"
		i__3 = jj - j + 1;
#line 493 "zlasyf.f"
		i__4 = *n - k;
#line 493 "zlasyf.f"
		z__1.r = -1., z__1.i = -0.;
#line 493 "zlasyf.f"
		zgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, 
			&a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 496 "zlasyf.f"
/* L40: */
#line 496 "zlasyf.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 500 "zlasyf.f"
	    i__2 = j - 1;
#line 500 "zlasyf.f"
	    i__3 = *n - k;
#line 500 "zlasyf.f"
	    z__1.r = -1., z__1.i = -0.;
#line 500 "zlasyf.f"
	    zgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, &a[(
		    k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * w_dim1], ldw,
		     &c_b1, &a[j * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
#line 503 "zlasyf.f"
/* L50: */
#line 503 "zlasyf.f"
	}

/*        Put U12 in standard form by partially undoing the interchanges */
/*        in columns k+1:n looping backwards from k+1 to n */

#line 508 "zlasyf.f"
	j = k + 1;
#line 509 "zlasyf.f"
L60:

/*           Undo the interchanges (if any) of rows JJ and JP at each */
/*           step J */

/*           (Here, J is a diagonal index) */
#line 515 "zlasyf.f"
	jj = j;
#line 516 "zlasyf.f"
	jp = ipiv[j];
#line 517 "zlasyf.f"
	if (jp < 0) {
#line 518 "zlasyf.f"
	    jp = -jp;
/*              (Here, J is a diagonal index) */
#line 520 "zlasyf.f"
	    ++j;
#line 521 "zlasyf.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length N-J+1 */
/*           of the rows to swap back doesn't include diagonal element) */
#line 524 "zlasyf.f"
	++j;
#line 525 "zlasyf.f"
	if (jp != jj && j <= *n) {
#line 525 "zlasyf.f"
	    i__1 = *n - j + 1;
#line 525 "zlasyf.f"
	    zswap_(&i__1, &a[jp + j * a_dim1], lda, &a[jj + j * a_dim1], lda);
#line 525 "zlasyf.f"
	}
#line 527 "zlasyf.f"
	if (j < *n) {
#line 527 "zlasyf.f"
	    goto L60;
#line 527 "zlasyf.f"
	}

/*        Set KB to the number of columns factorized */

#line 532 "zlasyf.f"
	*kb = *n - k;

#line 534 "zlasyf.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 */

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 542 "zlasyf.f"
	k = 1;
#line 543 "zlasyf.f"
L70:

/*        Exit from loop */

#line 547 "zlasyf.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 547 "zlasyf.f"
	    goto L90;
#line 547 "zlasyf.f"
	}

/*        Copy column K of A to column K of W and update it */

#line 552 "zlasyf.f"
	i__1 = *n - k + 1;
#line 552 "zlasyf.f"
	zcopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
#line 553 "zlasyf.f"
	i__1 = *n - k + 1;
#line 553 "zlasyf.f"
	i__2 = k - 1;
#line 553 "zlasyf.f"
	z__1.r = -1., z__1.i = -0.;
#line 553 "zlasyf.f"
	zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[k 
		+ w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1, (ftnlen)12);

#line 556 "zlasyf.f"
	kstep = 1;

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 561 "zlasyf.f"
	i__1 = k + k * w_dim1;
#line 561 "zlasyf.f"
	absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[k + k * 
		w_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */

#line 566 "zlasyf.f"
	if (k < *n) {
#line 567 "zlasyf.f"
	    i__1 = *n - k;
#line 567 "zlasyf.f"
	    imax = k + izamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 568 "zlasyf.f"
	    i__1 = imax + k * w_dim1;
#line 568 "zlasyf.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    k * w_dim1]), abs(d__2));
#line 569 "zlasyf.f"
	} else {
#line 570 "zlasyf.f"
	    colmax = 0.;
#line 571 "zlasyf.f"
	}

#line 573 "zlasyf.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 577 "zlasyf.f"
	    if (*info == 0) {
#line 577 "zlasyf.f"
		*info = k;
#line 577 "zlasyf.f"
	    }
#line 579 "zlasyf.f"
	    kp = k;
#line 580 "zlasyf.f"
	} else {
#line 581 "zlasyf.f"
	    if (absakk >= alpha * colmax) {

/*              no interchange, use 1-by-1 pivot block */

#line 585 "zlasyf.f"
		kp = k;
#line 586 "zlasyf.f"
	    } else {

/*              Copy column IMAX to column K+1 of W and update it */

#line 590 "zlasyf.f"
		i__1 = imax - k;
#line 590 "zlasyf.f"
		zcopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 591 "zlasyf.f"
		i__1 = *n - imax + 1;
#line 591 "zlasyf.f"
		zcopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 
			1) * w_dim1], &c__1);
#line 593 "zlasyf.f"
		i__1 = *n - k + 1;
#line 593 "zlasyf.f"
		i__2 = k - 1;
#line 593 "zlasyf.f"
		z__1.r = -1., z__1.i = -0.;
#line 593 "zlasyf.f"
		zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], 
			lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 1) * 
			w_dim1], &c__1, (ftnlen)12);

/*              JMAX is the column-index of the largest off-diagonal */
/*              element in row IMAX, and ROWMAX is its absolute value */

#line 600 "zlasyf.f"
		i__1 = imax - k;
#line 600 "zlasyf.f"
		jmax = k - 1 + izamax_(&i__1, &w[k + (k + 1) * w_dim1], &c__1)
			;
#line 601 "zlasyf.f"
		i__1 = jmax + (k + 1) * w_dim1;
#line 601 "zlasyf.f"
		rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			jmax + (k + 1) * w_dim1]), abs(d__2));
#line 602 "zlasyf.f"
		if (imax < *n) {
#line 603 "zlasyf.f"
		    i__1 = *n - imax;
#line 603 "zlasyf.f"
		    jmax = imax + izamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
/* Computing MAX */
#line 604 "zlasyf.f"
		    i__1 = jmax + (k + 1) * w_dim1;
#line 604 "zlasyf.f"
		    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) + (
			    d__2 = d_imag(&w[jmax + (k + 1) * w_dim1]), abs(
			    d__2));
#line 604 "zlasyf.f"
		    rowmax = max(d__3,d__4);
#line 605 "zlasyf.f"
		}

#line 607 "zlasyf.f"
		if (absakk >= alpha * colmax * (colmax / rowmax)) {

/*                 no interchange, use 1-by-1 pivot block */

#line 611 "zlasyf.f"
		    kp = k;
#line 612 "zlasyf.f"
		} else /* if(complicated condition) */ {
#line 612 "zlasyf.f"
		    i__1 = imax + (k + 1) * w_dim1;
#line 612 "zlasyf.f"
		    if ((d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    imax + (k + 1) * w_dim1]), abs(d__2)) >= alpha * 
			    rowmax) {

/*                 interchange rows and columns K and IMAX, use 1-by-1 */
/*                 pivot block */

#line 617 "zlasyf.f"
			kp = imax;

/*                 copy column K+1 of W to column K of W */

#line 621 "zlasyf.f"
			i__1 = *n - k + 1;
#line 621 "zlasyf.f"
			zcopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + 
				k * w_dim1], &c__1);
#line 622 "zlasyf.f"
		    } else {

/*                 interchange rows and columns K+1 and IMAX, use 2-by-2 */
/*                 pivot block */

#line 627 "zlasyf.f"
			kp = imax;
#line 628 "zlasyf.f"
			kstep = 2;
#line 629 "zlasyf.f"
		    }
#line 629 "zlasyf.f"
		}
#line 630 "zlasyf.f"
	    }

/*           ============================================================ */

/*           KK is the column of A where pivoting step stopped */

#line 636 "zlasyf.f"
	    kk = k + kstep - 1;

/*           Interchange rows and columns KP and KK. */
/*           Updated column KP is already stored in column KK of W. */

#line 641 "zlasyf.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP of submatrix A */
/*              at step K. No need to copy element into column K */
/*              (or K and K+1 for 2-by-2 pivot) of A, since these columns */
/*              will be later overwritten. */

#line 648 "zlasyf.f"
		i__1 = kp + kp * a_dim1;
#line 648 "zlasyf.f"
		i__2 = kk + kk * a_dim1;
#line 648 "zlasyf.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 649 "zlasyf.f"
		i__1 = kp - kk - 1;
#line 649 "zlasyf.f"
		zcopy_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 
			1) * a_dim1], lda);
#line 651 "zlasyf.f"
		if (kp < *n) {
#line 651 "zlasyf.f"
		    i__1 = *n - kp;
#line 651 "zlasyf.f"
		    zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 
			    + kp * a_dim1], &c__1);
#line 651 "zlasyf.f"
		}

/*              Interchange rows KK and KP in first K-1 columns of A */
/*              (columns K (or K and K+1 for 2-by-2 pivot) of A will be */
/*              later overwritten). Interchange rows KK and KP */
/*              in first KK columns of W. */

#line 659 "zlasyf.f"
		if (k > 1) {
#line 659 "zlasyf.f"
		    i__1 = k - 1;
#line 659 "zlasyf.f"
		    zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 659 "zlasyf.f"
		}
#line 661 "zlasyf.f"
		zswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 662 "zlasyf.f"
	    }

#line 664 "zlasyf.f"
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

#line 679 "zlasyf.f"
		i__1 = *n - k + 1;
#line 679 "zlasyf.f"
		zcopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 680 "zlasyf.f"
		if (k < *n) {
#line 681 "zlasyf.f"
		    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 681 "zlasyf.f"
		    r1.r = z__1.r, r1.i = z__1.i;
#line 682 "zlasyf.f"
		    i__1 = *n - k;
#line 682 "zlasyf.f"
		    zscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 683 "zlasyf.f"
		}

#line 685 "zlasyf.f"
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

#line 702 "zlasyf.f"
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

#line 729 "zlasyf.f"
		    i__1 = k + 1 + k * w_dim1;
#line 729 "zlasyf.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 730 "zlasyf.f"
		    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
#line 730 "zlasyf.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 731 "zlasyf.f"
		    z_div(&z__1, &w[k + k * w_dim1], &d21);
#line 731 "zlasyf.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 732 "zlasyf.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 732 "zlasyf.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 732 "zlasyf.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 732 "zlasyf.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 733 "zlasyf.f"
		    z_div(&z__1, &t, &d21);
#line 733 "zlasyf.f"
		    d21.r = z__1.r, d21.i = z__1.i;

/*                 Update elements in columns A(k) and A(k+1) as */
/*                 dot products of rows of ( W(k) W(k+1) ) and columns */
/*                 of D**(-1) */

#line 739 "zlasyf.f"
		    i__1 = *n;
#line 739 "zlasyf.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 740 "zlasyf.f"
			i__2 = j + k * a_dim1;
#line 740 "zlasyf.f"
			i__3 = j + k * w_dim1;
#line 740 "zlasyf.f"
			z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__3.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 740 "zlasyf.f"
			i__4 = j + (k + 1) * w_dim1;
#line 740 "zlasyf.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 740 "zlasyf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 740 "zlasyf.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 741 "zlasyf.f"
			i__2 = j + (k + 1) * a_dim1;
#line 741 "zlasyf.f"
			i__3 = j + (k + 1) * w_dim1;
#line 741 "zlasyf.f"
			z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__3.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 741 "zlasyf.f"
			i__4 = j + k * w_dim1;
#line 741 "zlasyf.f"
			z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4]
				.i;
#line 741 "zlasyf.f"
			z__1.r = d21.r * z__2.r - d21.i * z__2.i, z__1.i = 
				d21.r * z__2.i + d21.i * z__2.r;
#line 741 "zlasyf.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 742 "zlasyf.f"
/* L80: */
#line 742 "zlasyf.f"
		    }
#line 743 "zlasyf.f"
		}

/*              Copy D(k) to A */

#line 747 "zlasyf.f"
		i__1 = k + k * a_dim1;
#line 747 "zlasyf.f"
		i__2 = k + k * w_dim1;
#line 747 "zlasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 748 "zlasyf.f"
		i__1 = k + 1 + k * a_dim1;
#line 748 "zlasyf.f"
		i__2 = k + 1 + k * w_dim1;
#line 748 "zlasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 749 "zlasyf.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 749 "zlasyf.f"
		i__2 = k + 1 + (k + 1) * w_dim1;
#line 749 "zlasyf.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;

#line 751 "zlasyf.f"
	    }

#line 753 "zlasyf.f"
	}

/*        Store details of the interchanges in IPIV */

#line 757 "zlasyf.f"
	if (kstep == 1) {
#line 758 "zlasyf.f"
	    ipiv[k] = kp;
#line 759 "zlasyf.f"
	} else {
#line 760 "zlasyf.f"
	    ipiv[k] = -kp;
#line 761 "zlasyf.f"
	    ipiv[k + 1] = -kp;
#line 762 "zlasyf.f"
	}

/*        Increase K and return to the start of the main loop */

#line 766 "zlasyf.f"
	k += kstep;
#line 767 "zlasyf.f"
	goto L70;

#line 769 "zlasyf.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T */

/*        computing blocks of NB columns at a time */

#line 777 "zlasyf.f"
	i__1 = *n;
#line 777 "zlasyf.f"
	i__2 = *nb;
#line 777 "zlasyf.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 778 "zlasyf.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 778 "zlasyf.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 782 "zlasyf.f"
	    i__3 = j + jb - 1;
#line 782 "zlasyf.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 783 "zlasyf.f"
		i__4 = j + jb - jj;
#line 783 "zlasyf.f"
		i__5 = k - 1;
#line 783 "zlasyf.f"
		z__1.r = -1., z__1.i = -0.;
#line 783 "zlasyf.f"
		zgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1]
			, &c__1, (ftnlen)12);
#line 786 "zlasyf.f"
/* L100: */
#line 786 "zlasyf.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 790 "zlasyf.f"
	    if (j + jb <= *n) {
#line 790 "zlasyf.f"
		i__3 = *n - j - jb + 1;
#line 790 "zlasyf.f"
		i__4 = k - 1;
#line 790 "zlasyf.f"
		z__1.r = -1., z__1.i = -0.;
#line 790 "zlasyf.f"
		zgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, 
			&a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 790 "zlasyf.f"
	    }
#line 794 "zlasyf.f"
/* L110: */
#line 794 "zlasyf.f"
	}

/*        Put L21 in standard form by partially undoing the interchanges */
/*        of rows in columns 1:k-1 looping backwards from k-1 to 1 */

#line 799 "zlasyf.f"
	j = k - 1;
#line 800 "zlasyf.f"
L120:

/*           Undo the interchanges (if any) of rows JJ and JP at each */
/*           step J */

/*           (Here, J is a diagonal index) */
#line 806 "zlasyf.f"
	jj = j;
#line 807 "zlasyf.f"
	jp = ipiv[j];
#line 808 "zlasyf.f"
	if (jp < 0) {
#line 809 "zlasyf.f"
	    jp = -jp;
/*              (Here, J is a diagonal index) */
#line 811 "zlasyf.f"
	    --j;
#line 812 "zlasyf.f"
	}
/*           (NOTE: Here, J is used to determine row length. Length J */
/*           of the rows to swap back doesn't include diagonal element) */
#line 815 "zlasyf.f"
	--j;
#line 816 "zlasyf.f"
	if (jp != jj && j >= 1) {
#line 816 "zlasyf.f"
	    zswap_(&j, &a[jp + a_dim1], lda, &a[jj + a_dim1], lda);
#line 816 "zlasyf.f"
	}
#line 818 "zlasyf.f"
	if (j > 1) {
#line 818 "zlasyf.f"
	    goto L120;
#line 818 "zlasyf.f"
	}

/*        Set KB to the number of columns factorized */

#line 823 "zlasyf.f"
	*kb = k - 1;

#line 825 "zlasyf.f"
    }
#line 826 "zlasyf.f"
    return 0;

/*     End of ZLASYF */

} /* zlasyf_ */

