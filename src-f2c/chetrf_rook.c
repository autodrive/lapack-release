#line 1 "chetrf_rook.f"
/* chetrf_rook.f -- translated by f2c (version 20100827).
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

#line 1 "chetrf_rook.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CHETRF_ROOK computes the factorization of a complex Hermitian indefinite matrix using the bound
ed Bunch-Kaufman ("rook") diagonal pivoting method (blocked algorithm, calling Level 3 BLAS). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrf_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrf_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrf_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */

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
/* > CHETRF_ROOK computes the factorization of a comlex Hermitian matrix A */
/* > using the bounded Bunch-Kaufman ("rook") diagonal pivoting method. */
/* > The form of the factorization is */
/* > */
/* >    A = U*D*U**T  or  A = L*D*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is Hermitian and block diagonal with */
/* > 1-by-1 and 2-by-2 diagonal blocks. */
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
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)). */
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

/* > \date June 2016 */

/* > \ingroup complexHEcomputational */

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
/* >  June 2016,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int chetrf_rook__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer j, k, kb, nb;
    extern /* Subroutine */ int chetf2_rook__(char *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, ftnlen);
    static integer iws;
    extern /* Subroutine */ int clahef_rook__(char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin, iinfo;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.6.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 250 "chetrf_rook.f"
    /* Parameter adjustments */
#line 250 "chetrf_rook.f"
    a_dim1 = *lda;
#line 250 "chetrf_rook.f"
    a_offset = 1 + a_dim1;
#line 250 "chetrf_rook.f"
    a -= a_offset;
#line 250 "chetrf_rook.f"
    --ipiv;
#line 250 "chetrf_rook.f"
    --work;
#line 250 "chetrf_rook.f"

#line 250 "chetrf_rook.f"
    /* Function Body */
#line 250 "chetrf_rook.f"
    *info = 0;
#line 251 "chetrf_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 252 "chetrf_rook.f"
    lquery = *lwork == -1;
#line 253 "chetrf_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 254 "chetrf_rook.f"
	*info = -1;
#line 255 "chetrf_rook.f"
    } else if (*n < 0) {
#line 256 "chetrf_rook.f"
	*info = -2;
#line 257 "chetrf_rook.f"
    } else if (*lda < max(1,*n)) {
#line 258 "chetrf_rook.f"
	*info = -4;
#line 259 "chetrf_rook.f"
    } else if (*lwork < 1 && ! lquery) {
#line 260 "chetrf_rook.f"
	*info = -7;
#line 261 "chetrf_rook.f"
    }

#line 263 "chetrf_rook.f"
    if (*info == 0) {

/*        Determine the block size */

#line 267 "chetrf_rook.f"
	nb = ilaenv_(&c__1, "CHETRF_ROOK", uplo, n, &c_n1, &c_n1, &c_n1, (
		ftnlen)11, (ftnlen)1);
/* Computing MAX */
#line 268 "chetrf_rook.f"
	i__1 = 1, i__2 = *n * nb;
#line 268 "chetrf_rook.f"
	lwkopt = max(i__1,i__2);
#line 269 "chetrf_rook.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 270 "chetrf_rook.f"
    }

#line 272 "chetrf_rook.f"
    if (*info != 0) {
#line 273 "chetrf_rook.f"
	i__1 = -(*info);
#line 273 "chetrf_rook.f"
	xerbla_("CHETRF_ROOK", &i__1, (ftnlen)11);
#line 274 "chetrf_rook.f"
	return 0;
#line 275 "chetrf_rook.f"
    } else if (lquery) {
#line 276 "chetrf_rook.f"
	return 0;
#line 277 "chetrf_rook.f"
    }

#line 279 "chetrf_rook.f"
    nbmin = 2;
#line 280 "chetrf_rook.f"
    ldwork = *n;
#line 281 "chetrf_rook.f"
    if (nb > 1 && nb < *n) {
#line 282 "chetrf_rook.f"
	iws = ldwork * nb;
#line 283 "chetrf_rook.f"
	if (*lwork < iws) {
/* Computing MAX */
#line 284 "chetrf_rook.f"
	    i__1 = *lwork / ldwork;
#line 284 "chetrf_rook.f"
	    nb = max(i__1,1);
/* Computing MAX */
#line 285 "chetrf_rook.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CHETRF_ROOK", uplo, n, &c_n1, &
		    c_n1, &c_n1, (ftnlen)11, (ftnlen)1);
#line 285 "chetrf_rook.f"
	    nbmin = max(i__1,i__2);
#line 287 "chetrf_rook.f"
	}
#line 288 "chetrf_rook.f"
    } else {
#line 289 "chetrf_rook.f"
	iws = 1;
#line 290 "chetrf_rook.f"
    }
#line 291 "chetrf_rook.f"
    if (nb < nbmin) {
#line 291 "chetrf_rook.f"
	nb = *n;
#line 291 "chetrf_rook.f"
    }

#line 294 "chetrf_rook.f"
    if (upper) {

/*        Factorize A as U*D*U**T using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        KB, where KB is the number of columns factorized by CLAHEF_ROOK; */
/*        KB is either NB or NB-1, or K for the last block */

#line 302 "chetrf_rook.f"
	k = *n;
#line 303 "chetrf_rook.f"
L10:

/*        If K < 1, exit from loop */

#line 307 "chetrf_rook.f"
	if (k < 1) {
#line 307 "chetrf_rook.f"
	    goto L40;
#line 307 "chetrf_rook.f"
	}

#line 310 "chetrf_rook.f"
	if (k > nb) {

/*           Factorize columns k-kb+1:k of A and use blocked code to */
/*           update columns 1:k-kb */

#line 315 "chetrf_rook.f"
	    clahef_rook__(uplo, &k, &nb, &kb, &a[a_offset], lda, &ipiv[1], &
		    work[1], &ldwork, &iinfo, (ftnlen)1);
#line 317 "chetrf_rook.f"
	} else {

/*           Use unblocked code to factorize columns 1:k of A */

#line 321 "chetrf_rook.f"
	    chetf2_rook__(uplo, &k, &a[a_offset], lda, &ipiv[1], &iinfo, (
		    ftnlen)1);
#line 322 "chetrf_rook.f"
	    kb = k;
#line 323 "chetrf_rook.f"
	}

/*        Set INFO on the first occurrence of a zero pivot */

#line 327 "chetrf_rook.f"
	if (*info == 0 && iinfo > 0) {
#line 327 "chetrf_rook.f"
	    *info = iinfo;
#line 327 "chetrf_rook.f"
	}

/*        No need to adjust IPIV */

/*        Decrease K and return to the start of the main loop */

#line 334 "chetrf_rook.f"
	k -= kb;
#line 335 "chetrf_rook.f"
	goto L10;

#line 337 "chetrf_rook.f"
    } else {

/*        Factorize A as L*D*L**T using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        KB, where KB is the number of columns factorized by CLAHEF_ROOK; */
/*        KB is either NB or NB-1, or N-K+1 for the last block */

#line 345 "chetrf_rook.f"
	k = 1;
#line 346 "chetrf_rook.f"
L20:

/*        If K > N, exit from loop */

#line 350 "chetrf_rook.f"
	if (k > *n) {
#line 350 "chetrf_rook.f"
	    goto L40;
#line 350 "chetrf_rook.f"
	}

#line 353 "chetrf_rook.f"
	if (k <= *n - nb) {

/*           Factorize columns k:k+kb-1 of A and use blocked code to */
/*           update columns k+kb:n */

#line 358 "chetrf_rook.f"
	    i__1 = *n - k + 1;
#line 358 "chetrf_rook.f"
	    clahef_rook__(uplo, &i__1, &nb, &kb, &a[k + k * a_dim1], lda, &
		    ipiv[k], &work[1], &ldwork, &iinfo, (ftnlen)1);
#line 360 "chetrf_rook.f"
	} else {

/*           Use unblocked code to factorize columns k:n of A */

#line 364 "chetrf_rook.f"
	    i__1 = *n - k + 1;
#line 364 "chetrf_rook.f"
	    chetf2_rook__(uplo, &i__1, &a[k + k * a_dim1], lda, &ipiv[k], &
		    iinfo, (ftnlen)1);
#line 366 "chetrf_rook.f"
	    kb = *n - k + 1;
#line 367 "chetrf_rook.f"
	}

/*        Set INFO on the first occurrence of a zero pivot */

#line 371 "chetrf_rook.f"
	if (*info == 0 && iinfo > 0) {
#line 371 "chetrf_rook.f"
	    *info = iinfo + k - 1;
#line 371 "chetrf_rook.f"
	}

/*        Adjust IPIV */

#line 376 "chetrf_rook.f"
	i__1 = k + kb - 1;
#line 376 "chetrf_rook.f"
	for (j = k; j <= i__1; ++j) {
#line 377 "chetrf_rook.f"
	    if (ipiv[j] > 0) {
#line 378 "chetrf_rook.f"
		ipiv[j] = ipiv[j] + k - 1;
#line 379 "chetrf_rook.f"
	    } else {
#line 380 "chetrf_rook.f"
		ipiv[j] = ipiv[j] - k + 1;
#line 381 "chetrf_rook.f"
	    }
#line 382 "chetrf_rook.f"
/* L30: */
#line 382 "chetrf_rook.f"
	}

/*        Increase K and return to the start of the main loop */

#line 386 "chetrf_rook.f"
	k += kb;
#line 387 "chetrf_rook.f"
	goto L20;

#line 389 "chetrf_rook.f"
    }

#line 391 "chetrf_rook.f"
L40:
#line 392 "chetrf_rook.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 393 "chetrf_rook.f"
    return 0;

/*     End of CHETRF_ROOK */

} /* chetrf_rook__ */

