#line 1 "csytrf.f"
/* csytrf.f -- translated by f2c (version 20100827).
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

#line 1 "csytrf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CSYTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTRF computes the factorization of a complex symmetric matrix A */
/* > using the Bunch-Kaufman diagonal pivoting method.  The form of the */
/* > factorization is */
/* > */
/* >    A = U*D*U**T  or  A = L*D*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is symmetric and block diagonal with */
/* > with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
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
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
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
/* >          If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >          interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* >          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and */
/* >          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* >          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) = */
/* >          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were */
/* >          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of WORK.  LWORK >=1.  For best performance */
/* >          LWORK >= N*NB, where NB is the block size returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization */
/* >                has been completed, but the block diagonal matrix D is */
/* >                exactly singular, and division by zero will occur if it */
/* >                is used to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
/* > */
/*  ===================================================================== */
/* Subroutine */ int csytrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer j, k, kb, nb, iws;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    static logical upper;
    extern /* Subroutine */ int csytf2_(char *, integer *, doublecomplex *, 
	    integer *, integer *, integer *, ftnlen), xerbla_(char *, integer 
	    *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int clasyf_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 220 "csytrf.f"
    /* Parameter adjustments */
#line 220 "csytrf.f"
    a_dim1 = *lda;
#line 220 "csytrf.f"
    a_offset = 1 + a_dim1;
#line 220 "csytrf.f"
    a -= a_offset;
#line 220 "csytrf.f"
    --ipiv;
#line 220 "csytrf.f"
    --work;
#line 220 "csytrf.f"

#line 220 "csytrf.f"
    /* Function Body */
#line 220 "csytrf.f"
    *info = 0;
#line 221 "csytrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 222 "csytrf.f"
    lquery = *lwork == -1;
#line 223 "csytrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 224 "csytrf.f"
	*info = -1;
#line 225 "csytrf.f"
    } else if (*n < 0) {
#line 226 "csytrf.f"
	*info = -2;
#line 227 "csytrf.f"
    } else if (*lda < max(1,*n)) {
#line 228 "csytrf.f"
	*info = -4;
#line 229 "csytrf.f"
    } else if (*lwork < 1 && ! lquery) {
#line 230 "csytrf.f"
	*info = -7;
#line 231 "csytrf.f"
    }

#line 233 "csytrf.f"
    if (*info == 0) {

/*        Determine the block size */

#line 237 "csytrf.f"
	nb = ilaenv_(&c__1, "CSYTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
#line 238 "csytrf.f"
	lwkopt = *n * nb;
#line 239 "csytrf.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 240 "csytrf.f"
    }

#line 242 "csytrf.f"
    if (*info != 0) {
#line 243 "csytrf.f"
	i__1 = -(*info);
#line 243 "csytrf.f"
	xerbla_("CSYTRF", &i__1, (ftnlen)6);
#line 244 "csytrf.f"
	return 0;
#line 245 "csytrf.f"
    } else if (lquery) {
#line 246 "csytrf.f"
	return 0;
#line 247 "csytrf.f"
    }

#line 249 "csytrf.f"
    nbmin = 2;
#line 250 "csytrf.f"
    ldwork = *n;
#line 251 "csytrf.f"
    if (nb > 1 && nb < *n) {
#line 252 "csytrf.f"
	iws = ldwork * nb;
#line 253 "csytrf.f"
	if (*lwork < iws) {
/* Computing MAX */
#line 254 "csytrf.f"
	    i__1 = *lwork / ldwork;
#line 254 "csytrf.f"
	    nb = max(i__1,1);
/* Computing MAX */
#line 255 "csytrf.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CSYTRF", uplo, n, &c_n1, &c_n1, &
		    c_n1, (ftnlen)6, (ftnlen)1);
#line 255 "csytrf.f"
	    nbmin = max(i__1,i__2);
#line 256 "csytrf.f"
	}
#line 257 "csytrf.f"
    } else {
#line 258 "csytrf.f"
	iws = 1;
#line 259 "csytrf.f"
    }
#line 260 "csytrf.f"
    if (nb < nbmin) {
#line 260 "csytrf.f"
	nb = *n;
#line 260 "csytrf.f"
    }

#line 263 "csytrf.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        KB, where KB is the number of columns factorized by CLASYF; */
/*        KB is either NB or NB-1, or K for the last block */

#line 271 "csytrf.f"
	k = *n;
#line 272 "csytrf.f"
L10:

/*        If K < 1, exit from loop */

#line 276 "csytrf.f"
	if (k < 1) {
#line 276 "csytrf.f"
	    goto L40;
#line 276 "csytrf.f"
	}

#line 279 "csytrf.f"
	if (k > nb) {

/*           Factorize columns k-kb+1:k of A and use blocked code to */
/*           update columns 1:k-kb */

#line 284 "csytrf.f"
	    clasyf_(uplo, &k, &nb, &kb, &a[a_offset], lda, &ipiv[1], &work[1],
		     n, &iinfo, (ftnlen)1);
#line 285 "csytrf.f"
	} else {

/*           Use unblocked code to factorize columns 1:k of A */

#line 289 "csytrf.f"
	    csytf2_(uplo, &k, &a[a_offset], lda, &ipiv[1], &iinfo, (ftnlen)1);
#line 290 "csytrf.f"
	    kb = k;
#line 291 "csytrf.f"
	}

/*        Set INFO on the first occurrence of a zero pivot */

#line 295 "csytrf.f"
	if (*info == 0 && iinfo > 0) {
#line 295 "csytrf.f"
	    *info = iinfo;
#line 295 "csytrf.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 300 "csytrf.f"
	k -= kb;
#line 301 "csytrf.f"
	goto L10;

#line 303 "csytrf.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        KB, where KB is the number of columns factorized by CLASYF; */
/*        KB is either NB or NB-1, or N-K+1 for the last block */

#line 311 "csytrf.f"
	k = 1;
#line 312 "csytrf.f"
L20:

/*        If K > N, exit from loop */

#line 316 "csytrf.f"
	if (k > *n) {
#line 316 "csytrf.f"
	    goto L40;
#line 316 "csytrf.f"
	}

#line 319 "csytrf.f"
	if (k <= *n - nb) {

/*           Factorize columns k:k+kb-1 of A and use blocked code to */
/*           update columns k+kb:n */

#line 324 "csytrf.f"
	    i__1 = *n - k + 1;
#line 324 "csytrf.f"
	    clasyf_(uplo, &i__1, &nb, &kb, &a[k + k * a_dim1], lda, &ipiv[k], 
		    &work[1], n, &iinfo, (ftnlen)1);
#line 326 "csytrf.f"
	} else {

/*           Use unblocked code to factorize columns k:n of A */

#line 330 "csytrf.f"
	    i__1 = *n - k + 1;
#line 330 "csytrf.f"
	    csytf2_(uplo, &i__1, &a[k + k * a_dim1], lda, &ipiv[k], &iinfo, (
		    ftnlen)1);
#line 331 "csytrf.f"
	    kb = *n - k + 1;
#line 332 "csytrf.f"
	}

/*        Set INFO on the first occurrence of a zero pivot */

#line 336 "csytrf.f"
	if (*info == 0 && iinfo > 0) {
#line 336 "csytrf.f"
	    *info = iinfo + k - 1;
#line 336 "csytrf.f"
	}

/*        Adjust IPIV */

#line 341 "csytrf.f"
	i__1 = k + kb - 1;
#line 341 "csytrf.f"
	for (j = k; j <= i__1; ++j) {
#line 342 "csytrf.f"
	    if (ipiv[j] > 0) {
#line 343 "csytrf.f"
		ipiv[j] = ipiv[j] + k - 1;
#line 344 "csytrf.f"
	    } else {
#line 345 "csytrf.f"
		ipiv[j] = ipiv[j] - k + 1;
#line 346 "csytrf.f"
	    }
#line 347 "csytrf.f"
/* L30: */
#line 347 "csytrf.f"
	}

/*        Increase K and return to the start of the main loop */

#line 351 "csytrf.f"
	k += kb;
#line 352 "csytrf.f"
	goto L20;

#line 354 "csytrf.f"
    }

#line 356 "csytrf.f"
L40:
#line 357 "csytrf.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 358 "csytrf.f"
    return 0;

/*     End of CSYTRF */

} /* csytrf_ */

