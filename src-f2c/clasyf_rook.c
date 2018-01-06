#line 1 "clasyf_rook.f"
/* clasyf_rook.f -- translated by f2c (version 20100827).
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

#line 1 "clasyf_rook.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CLASYF_ROOK computes a partial factorization of a complex symmetric matrix using the bounded Bu
nch-Kaufman ("rook") diagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLASYF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasyf_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasyf_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasyf_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLASYF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */

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
/* > CLASYF_ROOK computes a partial factorization of a complex symmetric */
/* > matrix A using the bounded Bunch-Kaufman ("rook") diagonal */
/* > pivoting method. The partial factorization has the form: */
/* > */
/* > A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or: */
/* >       ( 0  U22 ) (  0   D  ) ( U12**T U22**T ) */
/* > */
/* > A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L' */
/* >       ( L21  I ) (  0  A22 ) (  0       I    ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > */
/* > CLASYF_ROOK is an auxiliary routine called by CSYTRF_ROOK. It uses */
/* > blocked code (calling Level 3 BLAS) to update the submatrix */
/* > A11 (if UPLO = 'U') or A22 (if UPLO = 'L'). */
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
/* >  November 2013,     Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int clasyf_rook__(char *uplo, integer *n, integer *nb, 
	integer *kb, doublecomplex *a, integer *lda, integer *ipiv, 
	doublecomplex *w, integer *ldw, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, p;
    static doublecomplex t, r1, d11, d12, d21, d22;
    static integer jb, ii, jj, kk, kp, kw, jp1, jp2, kkw;
    static logical done;
    static integer imax, jmax;
    static doublereal alpha;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal sfmin;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer itemp;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer kstep;
    static doublereal stemp, absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
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

#line 239 "clasyf_rook.f"
    /* Parameter adjustments */
#line 239 "clasyf_rook.f"
    a_dim1 = *lda;
#line 239 "clasyf_rook.f"
    a_offset = 1 + a_dim1;
#line 239 "clasyf_rook.f"
    a -= a_offset;
#line 239 "clasyf_rook.f"
    --ipiv;
#line 239 "clasyf_rook.f"
    w_dim1 = *ldw;
#line 239 "clasyf_rook.f"
    w_offset = 1 + w_dim1;
#line 239 "clasyf_rook.f"
    w -= w_offset;
#line 239 "clasyf_rook.f"

#line 239 "clasyf_rook.f"
    /* Function Body */
#line 239 "clasyf_rook.f"
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

#line 243 "clasyf_rook.f"
    alpha = (sqrt(17.) + 1.) / 8.;

/*     Compute machine safe minimum */

#line 247 "clasyf_rook.f"
    sfmin = slamch_("S", (ftnlen)1);

#line 249 "clasyf_rook.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 */

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

#line 257 "clasyf_rook.f"
	k = *n;
#line 258 "clasyf_rook.f"
L10:

/*        KW is the column of W which corresponds to column K of A */

#line 262 "clasyf_rook.f"
	kw = *nb + k - *n;

/*        Exit from loop */

#line 266 "clasyf_rook.f"
	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
#line 266 "clasyf_rook.f"
	    goto L30;
#line 266 "clasyf_rook.f"
	}

#line 269 "clasyf_rook.f"
	kstep = 1;
#line 270 "clasyf_rook.f"
	p = k;

/*        Copy column K of A to column KW of W and update it */

#line 274 "clasyf_rook.f"
	ccopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
#line 275 "clasyf_rook.f"
	if (k < *n) {
#line 275 "clasyf_rook.f"
	    i__1 = *n - k;
#line 275 "clasyf_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 275 "clasyf_rook.f"
	    cgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * 
		    w_dim1 + 1], &c__1, (ftnlen)12);
#line 275 "clasyf_rook.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 282 "clasyf_rook.f"
	i__1 = k + kw * w_dim1;
#line 282 "clasyf_rook.f"
	absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[k + kw * 
		w_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 288 "clasyf_rook.f"
	if (k > 1) {
#line 289 "clasyf_rook.f"
	    i__1 = k - 1;
#line 289 "clasyf_rook.f"
	    imax = icamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
#line 290 "clasyf_rook.f"
	    i__1 = imax + kw * w_dim1;
#line 290 "clasyf_rook.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    kw * w_dim1]), abs(d__2));
#line 291 "clasyf_rook.f"
	} else {
#line 292 "clasyf_rook.f"
	    colmax = 0.;
#line 293 "clasyf_rook.f"
	}

#line 295 "clasyf_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 299 "clasyf_rook.f"
	    if (*info == 0) {
#line 299 "clasyf_rook.f"
		*info = k;
#line 299 "clasyf_rook.f"
	    }
#line 301 "clasyf_rook.f"
	    kp = k;
#line 302 "clasyf_rook.f"
	    ccopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
#line 303 "clasyf_rook.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 312 "clasyf_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 316 "clasyf_rook.f"
		kp = k;

#line 318 "clasyf_rook.f"
	    } else {

#line 320 "clasyf_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 324 "clasyf_rook.f"
L12:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column KW-1 of W and update it */

#line 331 "clasyf_rook.f"
		ccopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
#line 332 "clasyf_rook.f"
		i__1 = k - imax;
#line 332 "clasyf_rook.f"
		ccopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);

#line 335 "clasyf_rook.f"
		if (k < *n) {
#line 335 "clasyf_rook.f"
		    i__1 = *n - k;
#line 335 "clasyf_rook.f"
		    z__1.r = -1., z__1.i = -0.;
#line 335 "clasyf_rook.f"
		    cgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1, (
			    ftnlen)12);
#line 335 "clasyf_rook.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 344 "clasyf_rook.f"
		if (imax != k) {
#line 345 "clasyf_rook.f"
		    i__1 = k - imax;
#line 345 "clasyf_rook.f"
		    jmax = imax + icamax_(&i__1, &w[imax + 1 + (kw - 1) * 
			    w_dim1], &c__1);
#line 347 "clasyf_rook.f"
		    i__1 = jmax + (kw - 1) * w_dim1;
#line 347 "clasyf_rook.f"
		    rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    w[jmax + (kw - 1) * w_dim1]), abs(d__2));
#line 348 "clasyf_rook.f"
		} else {
#line 349 "clasyf_rook.f"
		    rowmax = 0.;
#line 350 "clasyf_rook.f"
		}

#line 352 "clasyf_rook.f"
		if (imax > 1) {
#line 353 "clasyf_rook.f"
		    i__1 = imax - 1;
#line 353 "clasyf_rook.f"
		    itemp = icamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
#line 354 "clasyf_rook.f"
		    i__1 = itemp + (kw - 1) * w_dim1;
#line 354 "clasyf_rook.f"
		    stemp = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    itemp + (kw - 1) * w_dim1]), abs(d__2));
#line 355 "clasyf_rook.f"
		    if (stemp > rowmax) {
#line 356 "clasyf_rook.f"
			rowmax = stemp;
#line 357 "clasyf_rook.f"
			jmax = itemp;
#line 358 "clasyf_rook.f"
		    }
#line 359 "clasyf_rook.f"
		}

/*                 Equivalent to testing for */
/*                 CABS1( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 365 "clasyf_rook.f"
		i__1 = imax + (kw - 1) * w_dim1;
#line 365 "clasyf_rook.f"
		if (! ((d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax 
			+ (kw - 1) * w_dim1]), abs(d__2)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 371 "clasyf_rook.f"
		    kp = imax;

/*                    copy column KW-1 of W to column KW of W */

#line 375 "clasyf_rook.f"
		    ccopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 377 "clasyf_rook.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 382 "clasyf_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 388 "clasyf_rook.f"
		    kp = imax;
#line 389 "clasyf_rook.f"
		    kstep = 2;
#line 390 "clasyf_rook.f"
		    done = TRUE_;
#line 391 "clasyf_rook.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 395 "clasyf_rook.f"
		    p = imax;
#line 396 "clasyf_rook.f"
		    colmax = rowmax;
#line 397 "clasyf_rook.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 401 "clasyf_rook.f"
		    ccopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

#line 403 "clasyf_rook.f"
		}

/*                 End pivot search loop body */

#line 407 "clasyf_rook.f"
		if (! done) {
#line 407 "clasyf_rook.f"
		    goto L12;
#line 407 "clasyf_rook.f"
		}

#line 409 "clasyf_rook.f"
	    }

/*           ============================================================ */

#line 413 "clasyf_rook.f"
	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

#line 417 "clasyf_rook.f"
	    kkw = *nb + kk - *n;

#line 419 "clasyf_rook.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

#line 423 "clasyf_rook.f"
		i__1 = k - p;
#line 423 "clasyf_rook.f"
		ccopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * 
			a_dim1], lda);
#line 424 "clasyf_rook.f"
		ccopy_(&p, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &
			c__1);

/*              Interchange rows K and P in last N-K+1 columns of A */
/*              and last N-K+2 columns of W */

#line 429 "clasyf_rook.f"
		i__1 = *n - k + 1;
#line 429 "clasyf_rook.f"
		cswap_(&i__1, &a[k + k * a_dim1], lda, &a[p + k * a_dim1], 
			lda);
#line 430 "clasyf_rook.f"
		i__1 = *n - kk + 1;
#line 430 "clasyf_rook.f"
		cswap_(&i__1, &w[k + kkw * w_dim1], ldw, &w[p + kkw * w_dim1],
			 ldw);
#line 431 "clasyf_rook.f"
	    }

/*           Updated column KP is already stored in column KKW of W */

#line 435 "clasyf_rook.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

#line 439 "clasyf_rook.f"
		i__1 = kp + k * a_dim1;
#line 439 "clasyf_rook.f"
		i__2 = kk + k * a_dim1;
#line 439 "clasyf_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 440 "clasyf_rook.f"
		i__1 = k - 1 - kp;
#line 440 "clasyf_rook.f"
		ccopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
#line 442 "clasyf_rook.f"
		ccopy_(&kp, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
			c__1);

/*              Interchange rows KK and KP in last N-KK+1 columns */
/*              of A and W */

#line 447 "clasyf_rook.f"
		i__1 = *n - kk + 1;
#line 447 "clasyf_rook.f"
		cswap_(&i__1, &a[kk + kk * a_dim1], lda, &a[kp + kk * a_dim1],
			 lda);
#line 448 "clasyf_rook.f"
		i__1 = *n - kk + 1;
#line 448 "clasyf_rook.f"
		cswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
#line 450 "clasyf_rook.f"
	    }

#line 452 "clasyf_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column KW of W now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Store U(k) in column k of A */

#line 462 "clasyf_rook.f"
		ccopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 463 "clasyf_rook.f"
		if (k > 1) {
#line 464 "clasyf_rook.f"
		    i__1 = k + k * a_dim1;
#line 464 "clasyf_rook.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + 
			    k * a_dim1]), abs(d__2)) >= sfmin) {
#line 465 "clasyf_rook.f"
			z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 465 "clasyf_rook.f"
			r1.r = z__1.r, r1.i = z__1.i;
#line 466 "clasyf_rook.f"
			i__1 = k - 1;
#line 466 "clasyf_rook.f"
			cscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
#line 467 "clasyf_rook.f"
		    } else /* if(complicated condition) */ {
#line 467 "clasyf_rook.f"
			i__1 = k + k * a_dim1;
#line 467 "clasyf_rook.f"
			if (a[i__1].r != 0. || a[i__1].i != 0.) {
#line 468 "clasyf_rook.f"
			    i__1 = k - 1;
#line 468 "clasyf_rook.f"
			    for (ii = 1; ii <= i__1; ++ii) {
#line 469 "clasyf_rook.f"
				i__2 = ii + k * a_dim1;
#line 469 "clasyf_rook.f"
				z_div(&z__1, &a[ii + k * a_dim1], &a[k + k * 
					a_dim1]);
#line 469 "clasyf_rook.f"
				a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 470 "clasyf_rook.f"
/* L14: */
#line 470 "clasyf_rook.f"
			    }
#line 471 "clasyf_rook.f"
			}
#line 471 "clasyf_rook.f"
		    }
#line 472 "clasyf_rook.f"
		}

#line 474 "clasyf_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns KW and KW-1 of W now */
/*              hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

#line 484 "clasyf_rook.f"
		if (k > 2) {

/*                 Store U(k) and U(k-1) in columns k and k-1 of A */

#line 488 "clasyf_rook.f"
		    i__1 = k - 1 + kw * w_dim1;
#line 488 "clasyf_rook.f"
		    d12.r = w[i__1].r, d12.i = w[i__1].i;
#line 489 "clasyf_rook.f"
		    z_div(&z__1, &w[k + kw * w_dim1], &d12);
#line 489 "clasyf_rook.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 490 "clasyf_rook.f"
		    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d12);
#line 490 "clasyf_rook.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 491 "clasyf_rook.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 491 "clasyf_rook.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 491 "clasyf_rook.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 491 "clasyf_rook.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 492 "clasyf_rook.f"
		    i__1 = k - 2;
#line 492 "clasyf_rook.f"
		    for (j = 1; j <= i__1; ++j) {
#line 493 "clasyf_rook.f"
			i__2 = j + (k - 1) * a_dim1;
#line 493 "clasyf_rook.f"
			i__3 = j + (kw - 1) * w_dim1;
#line 493 "clasyf_rook.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 493 "clasyf_rook.f"
			i__4 = j + kw * w_dim1;
#line 493 "clasyf_rook.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 493 "clasyf_rook.f"
			z_div(&z__2, &z__3, &d12);
#line 493 "clasyf_rook.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 493 "clasyf_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 495 "clasyf_rook.f"
			i__2 = j + k * a_dim1;
#line 495 "clasyf_rook.f"
			i__3 = j + kw * w_dim1;
#line 495 "clasyf_rook.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 495 "clasyf_rook.f"
			i__4 = j + (kw - 1) * w_dim1;
#line 495 "clasyf_rook.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 495 "clasyf_rook.f"
			z_div(&z__2, &z__3, &d12);
#line 495 "clasyf_rook.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 495 "clasyf_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 497 "clasyf_rook.f"
/* L20: */
#line 497 "clasyf_rook.f"
		    }
#line 498 "clasyf_rook.f"
		}

/*              Copy D(k) to A */

#line 502 "clasyf_rook.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 502 "clasyf_rook.f"
		i__2 = k - 1 + (kw - 1) * w_dim1;
#line 502 "clasyf_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 503 "clasyf_rook.f"
		i__1 = k - 1 + k * a_dim1;
#line 503 "clasyf_rook.f"
		i__2 = k - 1 + kw * w_dim1;
#line 503 "clasyf_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 504 "clasyf_rook.f"
		i__1 = k + k * a_dim1;
#line 504 "clasyf_rook.f"
		i__2 = k + kw * w_dim1;
#line 504 "clasyf_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 505 "clasyf_rook.f"
	    }
#line 506 "clasyf_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 510 "clasyf_rook.f"
	if (kstep == 1) {
#line 511 "clasyf_rook.f"
	    ipiv[k] = kp;
#line 512 "clasyf_rook.f"
	} else {
#line 513 "clasyf_rook.f"
	    ipiv[k] = -p;
#line 514 "clasyf_rook.f"
	    ipiv[k - 1] = -kp;
#line 515 "clasyf_rook.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 519 "clasyf_rook.f"
	k -= kstep;
#line 520 "clasyf_rook.f"
	goto L10;

#line 522 "clasyf_rook.f"
L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T */

/*        computing blocks of NB columns at a time */

#line 530 "clasyf_rook.f"
	i__1 = -(*nb);
#line 530 "clasyf_rook.f"
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
#line 531 "clasyf_rook.f"
	    i__2 = *nb, i__3 = k - j + 1;
#line 531 "clasyf_rook.f"
	    jb = min(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

#line 535 "clasyf_rook.f"
	    i__2 = j + jb - 1;
#line 535 "clasyf_rook.f"
	    for (jj = j; jj <= i__2; ++jj) {
#line 536 "clasyf_rook.f"
		i__3 = jj - j + 1;
#line 536 "clasyf_rook.f"
		i__4 = *n - k;
#line 536 "clasyf_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 536 "clasyf_rook.f"
		cgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, 
			&a[j + jj * a_dim1], &c__1, (ftnlen)12);
#line 539 "clasyf_rook.f"
/* L40: */
#line 539 "clasyf_rook.f"
	    }

/*           Update the rectangular superdiagonal block */

#line 543 "clasyf_rook.f"
	    if (j >= 2) {
#line 543 "clasyf_rook.f"
		i__2 = j - 1;
#line 543 "clasyf_rook.f"
		i__3 = *n - k;
#line 543 "clasyf_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 543 "clasyf_rook.f"
		cgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, 
			&a[(k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * 
			w_dim1], ldw, &c_b1, &a[j * a_dim1 + 1], lda, (ftnlen)
			12, (ftnlen)9);
#line 543 "clasyf_rook.f"
	    }
#line 547 "clasyf_rook.f"
/* L50: */
#line 547 "clasyf_rook.f"
	}

/*        Put U12 in standard form by partially undoing the interchanges */
/*        in columns k+1:n */

#line 552 "clasyf_rook.f"
	j = k + 1;
#line 553 "clasyf_rook.f"
L60:

#line 555 "clasyf_rook.f"
	kstep = 1;
#line 556 "clasyf_rook.f"
	jp1 = 1;
#line 557 "clasyf_rook.f"
	jj = j;
#line 558 "clasyf_rook.f"
	jp2 = ipiv[j];
#line 559 "clasyf_rook.f"
	if (jp2 < 0) {
#line 560 "clasyf_rook.f"
	    jp2 = -jp2;
#line 561 "clasyf_rook.f"
	    ++j;
#line 562 "clasyf_rook.f"
	    jp1 = -ipiv[j];
#line 563 "clasyf_rook.f"
	    kstep = 2;
#line 564 "clasyf_rook.f"
	}

#line 566 "clasyf_rook.f"
	++j;
#line 567 "clasyf_rook.f"
	if (jp2 != jj && j <= *n) {
#line 567 "clasyf_rook.f"
	    i__1 = *n - j + 1;
#line 567 "clasyf_rook.f"
	    cswap_(&i__1, &a[jp2 + j * a_dim1], lda, &a[jj + j * a_dim1], lda)
		    ;
#line 567 "clasyf_rook.f"
	}
#line 569 "clasyf_rook.f"
	jj = j - 1;
#line 570 "clasyf_rook.f"
	if (jp1 != jj && kstep == 2) {
#line 570 "clasyf_rook.f"
	    i__1 = *n - j + 1;
#line 570 "clasyf_rook.f"
	    cswap_(&i__1, &a[jp1 + j * a_dim1], lda, &a[jj + j * a_dim1], lda)
		    ;
#line 570 "clasyf_rook.f"
	}
#line 572 "clasyf_rook.f"
	if (j <= *n) {
#line 572 "clasyf_rook.f"
	    goto L60;
#line 572 "clasyf_rook.f"
	}

/*        Set KB to the number of columns factorized */

#line 577 "clasyf_rook.f"
	*kb = *n - k;

#line 579 "clasyf_rook.f"
    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 */

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

#line 587 "clasyf_rook.f"
	k = 1;
#line 588 "clasyf_rook.f"
L70:

/*        Exit from loop */

#line 592 "clasyf_rook.f"
	if (k >= *nb && *nb < *n || k > *n) {
#line 592 "clasyf_rook.f"
	    goto L90;
#line 592 "clasyf_rook.f"
	}

#line 595 "clasyf_rook.f"
	kstep = 1;
#line 596 "clasyf_rook.f"
	p = k;

/*        Copy column K of A to column K of W and update it */

#line 600 "clasyf_rook.f"
	i__1 = *n - k + 1;
#line 600 "clasyf_rook.f"
	ccopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
#line 601 "clasyf_rook.f"
	if (k > 1) {
#line 601 "clasyf_rook.f"
	    i__1 = *n - k + 1;
#line 601 "clasyf_rook.f"
	    i__2 = k - 1;
#line 601 "clasyf_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 601 "clasyf_rook.f"
	    cgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &
		    w[k + w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1, (
		    ftnlen)12);
#line 601 "clasyf_rook.f"
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

#line 608 "clasyf_rook.f"
	i__1 = k + k * w_dim1;
#line 608 "clasyf_rook.f"
	absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[k + k * 
		w_dim1]), abs(d__2));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

#line 614 "clasyf_rook.f"
	if (k < *n) {
#line 615 "clasyf_rook.f"
	    i__1 = *n - k;
#line 615 "clasyf_rook.f"
	    imax = k + icamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
#line 616 "clasyf_rook.f"
	    i__1 = imax + k * w_dim1;
#line 616 "clasyf_rook.f"
	    colmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax + 
		    k * w_dim1]), abs(d__2));
#line 617 "clasyf_rook.f"
	} else {
#line 618 "clasyf_rook.f"
	    colmax = 0.;
#line 619 "clasyf_rook.f"
	}

#line 621 "clasyf_rook.f"
	if (max(absakk,colmax) == 0.) {

/*           Column K is zero or underflow: set INFO and continue */

#line 625 "clasyf_rook.f"
	    if (*info == 0) {
#line 625 "clasyf_rook.f"
		*info = k;
#line 625 "clasyf_rook.f"
	    }
#line 627 "clasyf_rook.f"
	    kp = k;
#line 628 "clasyf_rook.f"
	    i__1 = *n - k + 1;
#line 628 "clasyf_rook.f"
	    ccopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
		    c__1);
#line 629 "clasyf_rook.f"
	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

#line 638 "clasyf_rook.f"
	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

#line 642 "clasyf_rook.f"
		kp = k;

#line 644 "clasyf_rook.f"
	    } else {

#line 646 "clasyf_rook.f"
		done = FALSE_;

/*              Loop until pivot found */

#line 650 "clasyf_rook.f"
L72:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column K+1 of W and update it */

#line 657 "clasyf_rook.f"
		i__1 = imax - k;
#line 657 "clasyf_rook.f"
		ccopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
#line 658 "clasyf_rook.f"
		i__1 = *n - imax + 1;
#line 658 "clasyf_rook.f"
		ccopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 
			1) * w_dim1], &c__1);
#line 660 "clasyf_rook.f"
		if (k > 1) {
#line 660 "clasyf_rook.f"
		    i__1 = *n - k + 1;
#line 660 "clasyf_rook.f"
		    i__2 = k - 1;
#line 660 "clasyf_rook.f"
		    z__1.r = -1., z__1.i = -0.;
#line 660 "clasyf_rook.f"
		    cgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1]
			    , lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 
			    1) * w_dim1], &c__1, (ftnlen)12);
#line 660 "clasyf_rook.f"
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

#line 669 "clasyf_rook.f"
		if (imax != k) {
#line 670 "clasyf_rook.f"
		    i__1 = imax - k;
#line 670 "clasyf_rook.f"
		    jmax = k - 1 + icamax_(&i__1, &w[k + (k + 1) * w_dim1], &
			    c__1);
#line 671 "clasyf_rook.f"
		    i__1 = jmax + (k + 1) * w_dim1;
#line 671 "clasyf_rook.f"
		    rowmax = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&
			    w[jmax + (k + 1) * w_dim1]), abs(d__2));
#line 672 "clasyf_rook.f"
		} else {
#line 673 "clasyf_rook.f"
		    rowmax = 0.;
#line 674 "clasyf_rook.f"
		}

#line 676 "clasyf_rook.f"
		if (imax < *n) {
#line 677 "clasyf_rook.f"
		    i__1 = *n - imax;
#line 677 "clasyf_rook.f"
		    itemp = imax + icamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
#line 678 "clasyf_rook.f"
		    i__1 = itemp + (k + 1) * w_dim1;
#line 678 "clasyf_rook.f"
		    stemp = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[
			    itemp + (k + 1) * w_dim1]), abs(d__2));
#line 679 "clasyf_rook.f"
		    if (stemp > rowmax) {
#line 680 "clasyf_rook.f"
			rowmax = stemp;
#line 681 "clasyf_rook.f"
			jmax = itemp;
#line 682 "clasyf_rook.f"
		    }
#line 683 "clasyf_rook.f"
		}

/*                 Equivalent to testing for */
/*                 CABS1( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

#line 689 "clasyf_rook.f"
		i__1 = imax + (k + 1) * w_dim1;
#line 689 "clasyf_rook.f"
		if (! ((d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_imag(&w[imax 
			+ (k + 1) * w_dim1]), abs(d__2)) < alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

#line 695 "clasyf_rook.f"
		    kp = imax;

/*                    copy column K+1 of W to column K of W */

#line 699 "clasyf_rook.f"
		    i__1 = *n - k + 1;
#line 699 "clasyf_rook.f"
		    ccopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 701 "clasyf_rook.f"
		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

#line 706 "clasyf_rook.f"
		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

#line 712 "clasyf_rook.f"
		    kp = imax;
#line 713 "clasyf_rook.f"
		    kstep = 2;
#line 714 "clasyf_rook.f"
		    done = TRUE_;
#line 715 "clasyf_rook.f"
		} else {

/*                    Pivot not found: set params and repeat */

#line 719 "clasyf_rook.f"
		    p = imax;
#line 720 "clasyf_rook.f"
		    colmax = rowmax;
#line 721 "clasyf_rook.f"
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

#line 725 "clasyf_rook.f"
		    i__1 = *n - k + 1;
#line 725 "clasyf_rook.f"
		    ccopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

#line 727 "clasyf_rook.f"
		}

/*                 End pivot search loop body */

#line 731 "clasyf_rook.f"
		if (! done) {
#line 731 "clasyf_rook.f"
		    goto L72;
#line 731 "clasyf_rook.f"
		}

#line 733 "clasyf_rook.f"
	    }

/*           ============================================================ */

#line 737 "clasyf_rook.f"
	    kk = k + kstep - 1;

#line 739 "clasyf_rook.f"
	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

#line 743 "clasyf_rook.f"
		i__1 = p - k;
#line 743 "clasyf_rook.f"
		ccopy_(&i__1, &a[k + k * a_dim1], &c__1, &a[p + k * a_dim1], 
			lda);
#line 744 "clasyf_rook.f"
		i__1 = *n - p + 1;
#line 744 "clasyf_rook.f"
		ccopy_(&i__1, &a[p + k * a_dim1], &c__1, &a[p + p * a_dim1], &
			c__1);

/*              Interchange rows K and P in first K columns of A */
/*              and first K+1 columns of W */

#line 749 "clasyf_rook.f"
		cswap_(&k, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
#line 750 "clasyf_rook.f"
		cswap_(&kk, &w[k + w_dim1], ldw, &w[p + w_dim1], ldw);
#line 751 "clasyf_rook.f"
	    }

/*           Updated column KP is already stored in column KK of W */

#line 755 "clasyf_rook.f"
	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

#line 759 "clasyf_rook.f"
		i__1 = kp + k * a_dim1;
#line 759 "clasyf_rook.f"
		i__2 = kk + k * a_dim1;
#line 759 "clasyf_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 760 "clasyf_rook.f"
		i__1 = kp - k - 1;
#line 760 "clasyf_rook.f"
		ccopy_(&i__1, &a[k + 1 + kk * a_dim1], &c__1, &a[kp + (k + 1) 
			* a_dim1], lda);
#line 761 "clasyf_rook.f"
		i__1 = *n - kp + 1;
#line 761 "clasyf_rook.f"
		ccopy_(&i__1, &a[kp + kk * a_dim1], &c__1, &a[kp + kp * 
			a_dim1], &c__1);

/*              Interchange rows KK and KP in first KK columns of A and W */

#line 765 "clasyf_rook.f"
		cswap_(&kk, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
#line 766 "clasyf_rook.f"
		cswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
#line 767 "clasyf_rook.f"
	    }

#line 769 "clasyf_rook.f"
	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of W now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

/*              Store L(k) in column k of A */

#line 779 "clasyf_rook.f"
		i__1 = *n - k + 1;
#line 779 "clasyf_rook.f"
		ccopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
#line 780 "clasyf_rook.f"
		if (k < *n) {
#line 781 "clasyf_rook.f"
		    i__1 = k + k * a_dim1;
#line 781 "clasyf_rook.f"
		    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + 
			    k * a_dim1]), abs(d__2)) >= sfmin) {
#line 782 "clasyf_rook.f"
			z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 782 "clasyf_rook.f"
			r1.r = z__1.r, r1.i = z__1.i;
#line 783 "clasyf_rook.f"
			i__1 = *n - k;
#line 783 "clasyf_rook.f"
			cscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
#line 784 "clasyf_rook.f"
		    } else /* if(complicated condition) */ {
#line 784 "clasyf_rook.f"
			i__1 = k + k * a_dim1;
#line 784 "clasyf_rook.f"
			if (a[i__1].r != 0. || a[i__1].i != 0.) {
#line 785 "clasyf_rook.f"
			    i__1 = *n;
#line 785 "clasyf_rook.f"
			    for (ii = k + 1; ii <= i__1; ++ii) {
#line 786 "clasyf_rook.f"
				i__2 = ii + k * a_dim1;
#line 786 "clasyf_rook.f"
				z_div(&z__1, &a[ii + k * a_dim1], &a[k + k * 
					a_dim1]);
#line 786 "clasyf_rook.f"
				a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 787 "clasyf_rook.f"
/* L74: */
#line 787 "clasyf_rook.f"
			    }
#line 788 "clasyf_rook.f"
			}
#line 788 "clasyf_rook.f"
		    }
#line 789 "clasyf_rook.f"
		}

#line 791 "clasyf_rook.f"
	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

#line 800 "clasyf_rook.f"
		if (k < *n - 1) {

/*                 Store L(k) and L(k+1) in columns k and k+1 of A */

#line 804 "clasyf_rook.f"
		    i__1 = k + 1 + k * w_dim1;
#line 804 "clasyf_rook.f"
		    d21.r = w[i__1].r, d21.i = w[i__1].i;
#line 805 "clasyf_rook.f"
		    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
#line 805 "clasyf_rook.f"
		    d11.r = z__1.r, d11.i = z__1.i;
#line 806 "clasyf_rook.f"
		    z_div(&z__1, &w[k + k * w_dim1], &d21);
#line 806 "clasyf_rook.f"
		    d22.r = z__1.r, d22.i = z__1.i;
#line 807 "clasyf_rook.f"
		    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * 
			    d22.i + d11.i * d22.r;
#line 807 "clasyf_rook.f"
		    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 807 "clasyf_rook.f"
		    z_div(&z__1, &c_b1, &z__2);
#line 807 "clasyf_rook.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 808 "clasyf_rook.f"
		    i__1 = *n;
#line 808 "clasyf_rook.f"
		    for (j = k + 2; j <= i__1; ++j) {
#line 809 "clasyf_rook.f"
			i__2 = j + k * a_dim1;
#line 809 "clasyf_rook.f"
			i__3 = j + k * w_dim1;
#line 809 "clasyf_rook.f"
			z__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i, 
				z__4.i = d11.r * w[i__3].i + d11.i * w[i__3]
				.r;
#line 809 "clasyf_rook.f"
			i__4 = j + (k + 1) * w_dim1;
#line 809 "clasyf_rook.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 809 "clasyf_rook.f"
			z_div(&z__2, &z__3, &d21);
#line 809 "clasyf_rook.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 809 "clasyf_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 811 "clasyf_rook.f"
			i__2 = j + (k + 1) * a_dim1;
#line 811 "clasyf_rook.f"
			i__3 = j + (k + 1) * w_dim1;
#line 811 "clasyf_rook.f"
			z__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i, 
				z__4.i = d22.r * w[i__3].i + d22.i * w[i__3]
				.r;
#line 811 "clasyf_rook.f"
			i__4 = j + k * w_dim1;
#line 811 "clasyf_rook.f"
			z__3.r = z__4.r - w[i__4].r, z__3.i = z__4.i - w[i__4]
				.i;
#line 811 "clasyf_rook.f"
			z_div(&z__2, &z__3, &d21);
#line 811 "clasyf_rook.f"
			z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * 
				z__2.i + t.i * z__2.r;
#line 811 "clasyf_rook.f"
			a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 813 "clasyf_rook.f"
/* L80: */
#line 813 "clasyf_rook.f"
		    }
#line 814 "clasyf_rook.f"
		}

/*              Copy D(k) to A */

#line 818 "clasyf_rook.f"
		i__1 = k + k * a_dim1;
#line 818 "clasyf_rook.f"
		i__2 = k + k * w_dim1;
#line 818 "clasyf_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 819 "clasyf_rook.f"
		i__1 = k + 1 + k * a_dim1;
#line 819 "clasyf_rook.f"
		i__2 = k + 1 + k * w_dim1;
#line 819 "clasyf_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 820 "clasyf_rook.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 820 "clasyf_rook.f"
		i__2 = k + 1 + (k + 1) * w_dim1;
#line 820 "clasyf_rook.f"
		a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
#line 821 "clasyf_rook.f"
	    }
#line 822 "clasyf_rook.f"
	}

/*        Store details of the interchanges in IPIV */

#line 826 "clasyf_rook.f"
	if (kstep == 1) {
#line 827 "clasyf_rook.f"
	    ipiv[k] = kp;
#line 828 "clasyf_rook.f"
	} else {
#line 829 "clasyf_rook.f"
	    ipiv[k] = -p;
#line 830 "clasyf_rook.f"
	    ipiv[k + 1] = -kp;
#line 831 "clasyf_rook.f"
	}

/*        Increase K and return to the start of the main loop */

#line 835 "clasyf_rook.f"
	k += kstep;
#line 836 "clasyf_rook.f"
	goto L70;

#line 838 "clasyf_rook.f"
L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T */

/*        computing blocks of NB columns at a time */

#line 846 "clasyf_rook.f"
	i__1 = *n;
#line 846 "clasyf_rook.f"
	i__2 = *nb;
#line 846 "clasyf_rook.f"
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 847 "clasyf_rook.f"
	    i__3 = *nb, i__4 = *n - j + 1;
#line 847 "clasyf_rook.f"
	    jb = min(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

#line 851 "clasyf_rook.f"
	    i__3 = j + jb - 1;
#line 851 "clasyf_rook.f"
	    for (jj = j; jj <= i__3; ++jj) {
#line 852 "clasyf_rook.f"
		i__4 = j + jb - jj;
#line 852 "clasyf_rook.f"
		i__5 = k - 1;
#line 852 "clasyf_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 852 "clasyf_rook.f"
		cgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1]
			, &c__1, (ftnlen)12);
#line 855 "clasyf_rook.f"
/* L100: */
#line 855 "clasyf_rook.f"
	    }

/*           Update the rectangular subdiagonal block */

#line 859 "clasyf_rook.f"
	    if (j + jb <= *n) {
#line 859 "clasyf_rook.f"
		i__3 = *n - j - jb + 1;
#line 859 "clasyf_rook.f"
		i__4 = k - 1;
#line 859 "clasyf_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 859 "clasyf_rook.f"
		cgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, 
			&a[j + jb + j * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 859 "clasyf_rook.f"
	    }
#line 863 "clasyf_rook.f"
/* L110: */
#line 863 "clasyf_rook.f"
	}

/*        Put L21 in standard form by partially undoing the interchanges */
/*        in columns 1:k-1 */

#line 868 "clasyf_rook.f"
	j = k - 1;
#line 869 "clasyf_rook.f"
L120:

#line 871 "clasyf_rook.f"
	kstep = 1;
#line 872 "clasyf_rook.f"
	jp1 = 1;
#line 873 "clasyf_rook.f"
	jj = j;
#line 874 "clasyf_rook.f"
	jp2 = ipiv[j];
#line 875 "clasyf_rook.f"
	if (jp2 < 0) {
#line 876 "clasyf_rook.f"
	    jp2 = -jp2;
#line 877 "clasyf_rook.f"
	    --j;
#line 878 "clasyf_rook.f"
	    jp1 = -ipiv[j];
#line 879 "clasyf_rook.f"
	    kstep = 2;
#line 880 "clasyf_rook.f"
	}

#line 882 "clasyf_rook.f"
	--j;
#line 883 "clasyf_rook.f"
	if (jp2 != jj && j >= 1) {
#line 883 "clasyf_rook.f"
	    cswap_(&j, &a[jp2 + a_dim1], lda, &a[jj + a_dim1], lda);
#line 883 "clasyf_rook.f"
	}
#line 885 "clasyf_rook.f"
	jj = j + 1;
#line 886 "clasyf_rook.f"
	if (jp1 != jj && kstep == 2) {
#line 886 "clasyf_rook.f"
	    cswap_(&j, &a[jp1 + a_dim1], lda, &a[jj + a_dim1], lda);
#line 886 "clasyf_rook.f"
	}
#line 888 "clasyf_rook.f"
	if (j >= 1) {
#line 888 "clasyf_rook.f"
	    goto L120;
#line 888 "clasyf_rook.f"
	}

/*        Set KB to the number of columns factorized */

#line 893 "clasyf_rook.f"
	*kb = k - 1;

#line 895 "clasyf_rook.f"
    }
#line 896 "clasyf_rook.f"
    return 0;

/*     End of CLASYF_ROOK */

} /* clasyf_rook__ */

