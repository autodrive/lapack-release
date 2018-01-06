#line 1 "chetrf.f"
/* chetrf.f -- translated by f2c (version 20100827).
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

#line 1 "chetrf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CHETRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */

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
/* > CHETRF computes the factorization of a complex Hermitian matrix A */
/* > using the Bunch-Kaufman diagonal pivoting method.  The form of the */
/* > factorization is */
/* > */
/* >    A = U*D*U**H  or  A = L*D*L**H */
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

/* > \date December 2016 */

/* > \ingroup complexHEcomputational */

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
/* > */
/*  ===================================================================== */
/* Subroutine */ int chetrf_(char *uplo, integer *n, doublecomplex *a, 
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
    extern /* Subroutine */ int chetf2_(char *, integer *, doublecomplex *, 
	    integer *, integer *, integer *, ftnlen), clahef_(char *, integer 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 215 "chetrf.f"
    /* Parameter adjustments */
#line 215 "chetrf.f"
    a_dim1 = *lda;
#line 215 "chetrf.f"
    a_offset = 1 + a_dim1;
#line 215 "chetrf.f"
    a -= a_offset;
#line 215 "chetrf.f"
    --ipiv;
#line 215 "chetrf.f"
    --work;
#line 215 "chetrf.f"

#line 215 "chetrf.f"
    /* Function Body */
#line 215 "chetrf.f"
    *info = 0;
#line 216 "chetrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 217 "chetrf.f"
    lquery = *lwork == -1;
#line 218 "chetrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 219 "chetrf.f"
	*info = -1;
#line 220 "chetrf.f"
    } else if (*n < 0) {
#line 221 "chetrf.f"
	*info = -2;
#line 222 "chetrf.f"
    } else if (*lda < max(1,*n)) {
#line 223 "chetrf.f"
	*info = -4;
#line 224 "chetrf.f"
    } else if (*lwork < 1 && ! lquery) {
#line 225 "chetrf.f"
	*info = -7;
#line 226 "chetrf.f"
    }

#line 228 "chetrf.f"
    if (*info == 0) {

/*        Determine the block size */

#line 232 "chetrf.f"
	nb = ilaenv_(&c__1, "CHETRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
		 (ftnlen)1);
#line 233 "chetrf.f"
	lwkopt = *n * nb;
#line 234 "chetrf.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 235 "chetrf.f"
    }

#line 237 "chetrf.f"
    if (*info != 0) {
#line 238 "chetrf.f"
	i__1 = -(*info);
#line 238 "chetrf.f"
	xerbla_("CHETRF", &i__1, (ftnlen)6);
#line 239 "chetrf.f"
	return 0;
#line 240 "chetrf.f"
    } else if (lquery) {
#line 241 "chetrf.f"
	return 0;
#line 242 "chetrf.f"
    }

#line 244 "chetrf.f"
    nbmin = 2;
#line 245 "chetrf.f"
    ldwork = *n;
#line 246 "chetrf.f"
    if (nb > 1 && nb < *n) {
#line 247 "chetrf.f"
	iws = ldwork * nb;
#line 248 "chetrf.f"
	if (*lwork < iws) {
/* Computing MAX */
#line 249 "chetrf.f"
	    i__1 = *lwork / ldwork;
#line 249 "chetrf.f"
	    nb = max(i__1,1);
/* Computing MAX */
#line 250 "chetrf.f"
	    i__1 = 2, i__2 = ilaenv_(&c__2, "CHETRF", uplo, n, &c_n1, &c_n1, &
		    c_n1, (ftnlen)6, (ftnlen)1);
#line 250 "chetrf.f"
	    nbmin = max(i__1,i__2);
#line 251 "chetrf.f"
	}
#line 252 "chetrf.f"
    } else {
#line 253 "chetrf.f"
	iws = 1;
#line 254 "chetrf.f"
    }
#line 255 "chetrf.f"
    if (nb < nbmin) {
#line 255 "chetrf.f"
	nb = *n;
#line 255 "chetrf.f"
    }

#line 258 "chetrf.f"
    if (upper) {

/*        Factorize A as U*D*U**H using the upper triangle of A */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        KB, where KB is the number of columns factorized by CLAHEF; */
/*        KB is either NB or NB-1, or K for the last block */

#line 266 "chetrf.f"
	k = *n;
#line 267 "chetrf.f"
L10:

/*        If K < 1, exit from loop */

#line 271 "chetrf.f"
	if (k < 1) {
#line 271 "chetrf.f"
	    goto L40;
#line 271 "chetrf.f"
	}

#line 274 "chetrf.f"
	if (k > nb) {

/*           Factorize columns k-kb+1:k of A and use blocked code to */
/*           update columns 1:k-kb */

#line 279 "chetrf.f"
	    clahef_(uplo, &k, &nb, &kb, &a[a_offset], lda, &ipiv[1], &work[1],
		     n, &iinfo, (ftnlen)1);
#line 280 "chetrf.f"
	} else {

/*           Use unblocked code to factorize columns 1:k of A */

#line 284 "chetrf.f"
	    chetf2_(uplo, &k, &a[a_offset], lda, &ipiv[1], &iinfo, (ftnlen)1);
#line 285 "chetrf.f"
	    kb = k;
#line 286 "chetrf.f"
	}

/*        Set INFO on the first occurrence of a zero pivot */

#line 290 "chetrf.f"
	if (*info == 0 && iinfo > 0) {
#line 290 "chetrf.f"
	    *info = iinfo;
#line 290 "chetrf.f"
	}

/*        Decrease K and return to the start of the main loop */

#line 295 "chetrf.f"
	k -= kb;
#line 296 "chetrf.f"
	goto L10;

#line 298 "chetrf.f"
    } else {

/*        Factorize A as L*D*L**H using the lower triangle of A */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        KB, where KB is the number of columns factorized by CLAHEF; */
/*        KB is either NB or NB-1, or N-K+1 for the last block */

#line 306 "chetrf.f"
	k = 1;
#line 307 "chetrf.f"
L20:

/*        If K > N, exit from loop */

#line 311 "chetrf.f"
	if (k > *n) {
#line 311 "chetrf.f"
	    goto L40;
#line 311 "chetrf.f"
	}

#line 314 "chetrf.f"
	if (k <= *n - nb) {

/*           Factorize columns k:k+kb-1 of A and use blocked code to */
/*           update columns k+kb:n */

#line 319 "chetrf.f"
	    i__1 = *n - k + 1;
#line 319 "chetrf.f"
	    clahef_(uplo, &i__1, &nb, &kb, &a[k + k * a_dim1], lda, &ipiv[k], 
		    &work[1], n, &iinfo, (ftnlen)1);
#line 321 "chetrf.f"
	} else {

/*           Use unblocked code to factorize columns k:n of A */

#line 325 "chetrf.f"
	    i__1 = *n - k + 1;
#line 325 "chetrf.f"
	    chetf2_(uplo, &i__1, &a[k + k * a_dim1], lda, &ipiv[k], &iinfo, (
		    ftnlen)1);
#line 326 "chetrf.f"
	    kb = *n - k + 1;
#line 327 "chetrf.f"
	}

/*        Set INFO on the first occurrence of a zero pivot */

#line 331 "chetrf.f"
	if (*info == 0 && iinfo > 0) {
#line 331 "chetrf.f"
	    *info = iinfo + k - 1;
#line 331 "chetrf.f"
	}

/*        Adjust IPIV */

#line 336 "chetrf.f"
	i__1 = k + kb - 1;
#line 336 "chetrf.f"
	for (j = k; j <= i__1; ++j) {
#line 337 "chetrf.f"
	    if (ipiv[j] > 0) {
#line 338 "chetrf.f"
		ipiv[j] = ipiv[j] + k - 1;
#line 339 "chetrf.f"
	    } else {
#line 340 "chetrf.f"
		ipiv[j] = ipiv[j] - k + 1;
#line 341 "chetrf.f"
	    }
#line 342 "chetrf.f"
/* L30: */
#line 342 "chetrf.f"
	}

/*        Increase K and return to the start of the main loop */

#line 346 "chetrf.f"
	k += kb;
#line 347 "chetrf.f"
	goto L20;

#line 349 "chetrf.f"
    }

#line 351 "chetrf.f"
L40:
#line 352 "chetrf.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 353 "chetrf.f"
    return 0;

/*     End of CHETRF */

} /* chetrf_ */

