#line 1 "csytrf_aa.f"
/* csytrf_aa.f -- translated by f2c (version 20100827).
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

#line 1 "csytrf_aa.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublecomplex c_b15 = {1.,0.};
static doublecomplex c_b19 = {-1.,-0.};

/* > \brief \b CSYTRF_AA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRF_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrf_
aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrf_
aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrf_
aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            N, LDA, LWORK, INFO */
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
/* > CSYTRF_AA computes the factorization of a complex symmetric matrix A */
/* > using the Aasen's algorithm.  The form of the factorization is */
/* > */
/* >    A = U*T*U**T  or  A = L*T*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is a complex symmetric tridiagonal matrix. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, the tridiagonal matrix is stored in the diagonals */
/* >          and the subdiagonals of A just below (or above) the diagonals, */
/* >          and L is stored below (or above) the subdiaonals, when UPLO */
/* >          is 'L' (or 'U'). */
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
/* >          On exit, it contains the details of the interchanges, i.e., */
/* >          the row and column k of A were interchanged with the */
/* >          row and column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of WORK. LWORK >= MAX(1,2*N). For optimum performance */
/* >          LWORK >= N*(1+NB), where NB is the optimal blocksize. */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csytrf_aa__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int clasyf_aa__(char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, ftnlen);
    static integer j1, k1, k2, j2, j3, jb, nb, mj, nj;
    static doublecomplex alpha;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    cswap_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*     .. Parameters .. */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Determine the block size */

#line 178 "csytrf_aa.f"
    /* Parameter adjustments */
#line 178 "csytrf_aa.f"
    a_dim1 = *lda;
#line 178 "csytrf_aa.f"
    a_offset = 1 + a_dim1;
#line 178 "csytrf_aa.f"
    a -= a_offset;
#line 178 "csytrf_aa.f"
    --ipiv;
#line 178 "csytrf_aa.f"
    --work;
#line 178 "csytrf_aa.f"

#line 178 "csytrf_aa.f"
    /* Function Body */
#line 178 "csytrf_aa.f"
    nb = ilaenv_(&c__1, "CSYTRF_AA", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)9, 
	    (ftnlen)1);

/*     Test the input parameters. */

#line 182 "csytrf_aa.f"
    *info = 0;
#line 183 "csytrf_aa.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 184 "csytrf_aa.f"
    lquery = *lwork == -1;
#line 185 "csytrf_aa.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 186 "csytrf_aa.f"
	*info = -1;
#line 187 "csytrf_aa.f"
    } else if (*n < 0) {
#line 188 "csytrf_aa.f"
	*info = -2;
#line 189 "csytrf_aa.f"
    } else if (*lda < max(1,*n)) {
#line 190 "csytrf_aa.f"
	*info = -4;
#line 191 "csytrf_aa.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 191 "csytrf_aa.f"
	i__1 = 1, i__2 = *n << 1;
#line 191 "csytrf_aa.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 192 "csytrf_aa.f"
	    *info = -7;
#line 193 "csytrf_aa.f"
	}
#line 193 "csytrf_aa.f"
    }

#line 195 "csytrf_aa.f"
    if (*info == 0) {
#line 196 "csytrf_aa.f"
	lwkopt = (nb + 1) * *n;
#line 197 "csytrf_aa.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 198 "csytrf_aa.f"
    }

#line 200 "csytrf_aa.f"
    if (*info != 0) {
#line 201 "csytrf_aa.f"
	i__1 = -(*info);
#line 201 "csytrf_aa.f"
	xerbla_("CSYTRF_AA", &i__1, (ftnlen)9);
#line 202 "csytrf_aa.f"
	return 0;
#line 203 "csytrf_aa.f"
    } else if (lquery) {
#line 204 "csytrf_aa.f"
	return 0;
#line 205 "csytrf_aa.f"
    }

/*     Quick return */

#line 209 "csytrf_aa.f"
    if (*n == 0) {
#line 210 "csytrf_aa.f"
	return 0;
#line 211 "csytrf_aa.f"
    }
#line 212 "csytrf_aa.f"
    ipiv[1] = 1;
#line 213 "csytrf_aa.f"
    if (*n == 1) {
#line 214 "csytrf_aa.f"
	return 0;
#line 215 "csytrf_aa.f"
    }

/*     Adjust block size based on the workspace size */

#line 219 "csytrf_aa.f"
    if (*lwork < (nb + 1) * *n) {
#line 220 "csytrf_aa.f"
	nb = (*lwork - *n) / *n;
#line 221 "csytrf_aa.f"
    }

#line 223 "csytrf_aa.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the upper triangle of A */
/*        ..................................................... */

/*        Copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N)) */

#line 231 "csytrf_aa.f"
	ccopy_(n, &a[a_dim1 + 1], lda, &work[1], &c__1);

/*        J is the main loop index, increasing from 1 to N in steps of */
/*        JB, where JB is the number of columns factorized by CLASYF; */
/*        JB is either NB, or N-J+1 for the last block */

#line 237 "csytrf_aa.f"
	j = 0;
#line 238 "csytrf_aa.f"
L10:
#line 239 "csytrf_aa.f"
	if (j >= *n) {
#line 239 "csytrf_aa.f"
	    goto L20;
#line 239 "csytrf_aa.f"
	}

/*        each step of the main loop */
/*         J is the last column of the previous panel */
/*         J1 is the first column of the current panel */
/*         K1 identifies if the previous column of the panel has been */
/*          explicitly stored, e.g., K1=1 for the first panel, and */
/*          K1=0 for the rest */

#line 249 "csytrf_aa.f"
	j1 = j + 1;
/* Computing MIN */
#line 250 "csytrf_aa.f"
	i__1 = *n - j1 + 1;
#line 250 "csytrf_aa.f"
	jb = min(i__1,nb);
#line 251 "csytrf_aa.f"
	k1 = max(1,j) - j;

/*        Panel factorization */

#line 255 "csytrf_aa.f"
	i__1 = 2 - k1;
#line 255 "csytrf_aa.f"
	i__2 = *n - j;
#line 255 "csytrf_aa.f"
	clasyf_aa__(uplo, &i__1, &i__2, &jb, &a[max(1,j) + (j + 1) * a_dim1], 
		lda, &ipiv[j + 1], &work[1], n, &work[*n * nb + 1], (ftnlen)1)
		;

/*        Ajust IPIV and apply it back (J-th step picks (J+1)-th pivot) */

/* Computing MIN */
#line 261 "csytrf_aa.f"
	i__2 = *n, i__3 = j + jb + 1;
#line 261 "csytrf_aa.f"
	i__1 = min(i__2,i__3);
#line 261 "csytrf_aa.f"
	for (j2 = j + 2; j2 <= i__1; ++j2) {
#line 262 "csytrf_aa.f"
	    ipiv[j2] += j;
#line 263 "csytrf_aa.f"
	    if (j2 != ipiv[j2] && j1 - k1 > 2) {
#line 264 "csytrf_aa.f"
		i__2 = j1 - k1 - 2;
#line 264 "csytrf_aa.f"
		cswap_(&i__2, &a[j2 * a_dim1 + 1], &c__1, &a[ipiv[j2] * 
			a_dim1 + 1], &c__1);
#line 266 "csytrf_aa.f"
	    }
#line 267 "csytrf_aa.f"
	}
#line 268 "csytrf_aa.f"
	j += jb;

/*        Trailing submatrix update, where */
/*         the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and */
/*         WORK stores the current block of the auxiriarly matrix H */

#line 274 "csytrf_aa.f"
	if (j < *n) {

/*           If first panel and JB=1 (NB=1), then nothing to do */

#line 278 "csytrf_aa.f"
	    if (j1 > 1 || jb > 1) {

/*              Merge rank-1 update with BLAS-3 update */

#line 282 "csytrf_aa.f"
		i__1 = j + (j + 1) * a_dim1;
#line 282 "csytrf_aa.f"
		alpha.r = a[i__1].r, alpha.i = a[i__1].i;
#line 283 "csytrf_aa.f"
		i__1 = j + (j + 1) * a_dim1;
#line 283 "csytrf_aa.f"
		a[i__1].r = 1., a[i__1].i = 0.;
#line 284 "csytrf_aa.f"
		i__1 = *n - j;
#line 284 "csytrf_aa.f"
		ccopy_(&i__1, &a[j - 1 + (j + 1) * a_dim1], lda, &work[j + 1 
			- j1 + 1 + jb * *n], &c__1);
#line 286 "csytrf_aa.f"
		i__1 = *n - j;
#line 286 "csytrf_aa.f"
		cscal_(&i__1, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);

/*              K1 identifies if the previous column of the panel has been */
/*               explicitly stored, e.g., K1=1 and K2= 0 for the first panel, */
/*               while K1=0 and K2=1 for the rest */

#line 292 "csytrf_aa.f"
		if (j1 > 1) {

/*                 Not first panel */

#line 296 "csytrf_aa.f"
		    k2 = 1;
#line 297 "csytrf_aa.f"
		} else {

/*                 First panel */

#line 301 "csytrf_aa.f"
		    k2 = 0;

/*                 First update skips the first column */

#line 305 "csytrf_aa.f"
		    --jb;
#line 306 "csytrf_aa.f"
		}

#line 308 "csytrf_aa.f"
		i__1 = *n;
#line 308 "csytrf_aa.f"
		i__2 = nb;
#line 308 "csytrf_aa.f"
		for (j2 = j + 1; i__2 < 0 ? j2 >= i__1 : j2 <= i__1; j2 += 
			i__2) {
/* Computing MIN */
#line 309 "csytrf_aa.f"
		    i__3 = nb, i__4 = *n - j2 + 1;
#line 309 "csytrf_aa.f"
		    nj = min(i__3,i__4);

/*                 Update (J2, J2) diagonal block with CGEMV */

#line 313 "csytrf_aa.f"
		    j3 = j2;
#line 314 "csytrf_aa.f"
		    for (mj = nj - 1; mj >= 1; --mj) {
#line 315 "csytrf_aa.f"
			i__3 = jb + 1;
#line 315 "csytrf_aa.f"
			cgemv_("No transpose", &mj, &i__3, &c_b19, &work[j3 - 
				j1 + 1 + k1 * *n], n, &a[j1 - k2 + j3 * 
				a_dim1], &c__1, &c_b15, &a[j3 + j3 * a_dim1], 
				lda, (ftnlen)12);
#line 319 "csytrf_aa.f"
			++j3;
#line 320 "csytrf_aa.f"
		    }

/*                 Update off-diagonal block of J2-th block row with CGEMM */

#line 324 "csytrf_aa.f"
		    i__3 = *n - j3 + 1;
#line 324 "csytrf_aa.f"
		    i__4 = jb + 1;
#line 324 "csytrf_aa.f"
		    cgemm_("Transpose", "Transpose", &nj, &i__3, &i__4, &
			    c_b19, &a[j1 - k2 + j2 * a_dim1], lda, &work[j3 - 
			    j1 + 1 + k1 * *n], n, &c_b15, &a[j2 + j3 * a_dim1]
			    , lda, (ftnlen)9, (ftnlen)9);
#line 329 "csytrf_aa.f"
		}

/*              Recover T( J, J+1 ) */

#line 333 "csytrf_aa.f"
		i__2 = j + (j + 1) * a_dim1;
#line 333 "csytrf_aa.f"
		a[i__2].r = alpha.r, a[i__2].i = alpha.i;
#line 334 "csytrf_aa.f"
	    }

/*           WORK(J+1, 1) stores H(J+1, 1) */

#line 338 "csytrf_aa.f"
	    i__2 = *n - j;
#line 338 "csytrf_aa.f"
	    ccopy_(&i__2, &a[j + 1 + (j + 1) * a_dim1], lda, &work[1], &c__1);
#line 339 "csytrf_aa.f"
	}
#line 340 "csytrf_aa.f"
	goto L10;
#line 341 "csytrf_aa.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

/*        copy first column A(1:N, 1) into H(1:N, 1) */
/*         (stored in WORK(1:N)) */

#line 350 "csytrf_aa.f"
	ccopy_(n, &a[a_dim1 + 1], &c__1, &work[1], &c__1);

/*        J is the main loop index, increasing from 1 to N in steps of */
/*        JB, where JB is the number of columns factorized by CLASYF; */
/*        JB is either NB, or N-J+1 for the last block */

#line 356 "csytrf_aa.f"
	j = 0;
#line 357 "csytrf_aa.f"
L11:
#line 358 "csytrf_aa.f"
	if (j >= *n) {
#line 358 "csytrf_aa.f"
	    goto L20;
#line 358 "csytrf_aa.f"
	}

/*        each step of the main loop */
/*         J is the last column of the previous panel */
/*         J1 is the first column of the current panel */
/*         K1 identifies if the previous column of the panel has been */
/*          explicitly stored, e.g., K1=1 for the first panel, and */
/*          K1=0 for the rest */

#line 368 "csytrf_aa.f"
	j1 = j + 1;
/* Computing MIN */
#line 369 "csytrf_aa.f"
	i__2 = *n - j1 + 1;
#line 369 "csytrf_aa.f"
	jb = min(i__2,nb);
#line 370 "csytrf_aa.f"
	k1 = max(1,j) - j;

/*        Panel factorization */

#line 374 "csytrf_aa.f"
	i__2 = 2 - k1;
#line 374 "csytrf_aa.f"
	i__1 = *n - j;
#line 374 "csytrf_aa.f"
	clasyf_aa__(uplo, &i__2, &i__1, &jb, &a[j + 1 + max(1,j) * a_dim1], 
		lda, &ipiv[j + 1], &work[1], n, &work[*n * nb + 1], (ftnlen)1)
		;

/*        Ajust IPIV and apply it back (J-th step picks (J+1)-th pivot) */

/* Computing MIN */
#line 380 "csytrf_aa.f"
	i__1 = *n, i__3 = j + jb + 1;
#line 380 "csytrf_aa.f"
	i__2 = min(i__1,i__3);
#line 380 "csytrf_aa.f"
	for (j2 = j + 2; j2 <= i__2; ++j2) {
#line 381 "csytrf_aa.f"
	    ipiv[j2] += j;
#line 382 "csytrf_aa.f"
	    if (j2 != ipiv[j2] && j1 - k1 > 2) {
#line 383 "csytrf_aa.f"
		i__1 = j1 - k1 - 2;
#line 383 "csytrf_aa.f"
		cswap_(&i__1, &a[j2 + a_dim1], lda, &a[ipiv[j2] + a_dim1], 
			lda);
#line 385 "csytrf_aa.f"
	    }
#line 386 "csytrf_aa.f"
	}
#line 387 "csytrf_aa.f"
	j += jb;

/*        Trailing submatrix update, where */
/*          A(J2+1, J1-1) stores L(J2+1, J1) and */
/*          WORK(J2+1, 1) stores H(J2+1, 1) */

#line 393 "csytrf_aa.f"
	if (j < *n) {

/*           if first panel and JB=1 (NB=1), then nothing to do */

#line 397 "csytrf_aa.f"
	    if (j1 > 1 || jb > 1) {

/*              Merge rank-1 update with BLAS-3 update */

#line 401 "csytrf_aa.f"
		i__2 = j + 1 + j * a_dim1;
#line 401 "csytrf_aa.f"
		alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 402 "csytrf_aa.f"
		i__2 = j + 1 + j * a_dim1;
#line 402 "csytrf_aa.f"
		a[i__2].r = 1., a[i__2].i = 0.;
#line 403 "csytrf_aa.f"
		i__2 = *n - j;
#line 403 "csytrf_aa.f"
		ccopy_(&i__2, &a[j + 1 + (j - 1) * a_dim1], &c__1, &work[j + 
			1 - j1 + 1 + jb * *n], &c__1);
#line 405 "csytrf_aa.f"
		i__2 = *n - j;
#line 405 "csytrf_aa.f"
		cscal_(&i__2, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);

/*              K1 identifies if the previous column of the panel has been */
/*               explicitly stored, e.g., K1=1 and K2= 0 for the first panel, */
/*               while K1=0 and K2=1 for the rest */

#line 411 "csytrf_aa.f"
		if (j1 > 1) {

/*                 Not first panel */

#line 415 "csytrf_aa.f"
		    k2 = 1;
#line 416 "csytrf_aa.f"
		} else {

/*                 First panel */

#line 420 "csytrf_aa.f"
		    k2 = 0;

/*                 First update skips the first column */

#line 424 "csytrf_aa.f"
		    --jb;
#line 425 "csytrf_aa.f"
		}

#line 427 "csytrf_aa.f"
		i__2 = *n;
#line 427 "csytrf_aa.f"
		i__1 = nb;
#line 427 "csytrf_aa.f"
		for (j2 = j + 1; i__1 < 0 ? j2 >= i__2 : j2 <= i__2; j2 += 
			i__1) {
/* Computing MIN */
#line 428 "csytrf_aa.f"
		    i__3 = nb, i__4 = *n - j2 + 1;
#line 428 "csytrf_aa.f"
		    nj = min(i__3,i__4);

/*                 Update (J2, J2) diagonal block with CGEMV */

#line 432 "csytrf_aa.f"
		    j3 = j2;
#line 433 "csytrf_aa.f"
		    for (mj = nj - 1; mj >= 1; --mj) {
#line 434 "csytrf_aa.f"
			i__3 = jb + 1;
#line 434 "csytrf_aa.f"
			cgemv_("No transpose", &mj, &i__3, &c_b19, &work[j3 - 
				j1 + 1 + k1 * *n], n, &a[j3 + (j1 - k2) * 
				a_dim1], lda, &c_b15, &a[j3 + j3 * a_dim1], &
				c__1, (ftnlen)12);
#line 438 "csytrf_aa.f"
			++j3;
#line 439 "csytrf_aa.f"
		    }

/*                 Update off-diagonal block in J2-th block column with CGEMM */

#line 443 "csytrf_aa.f"
		    i__3 = *n - j3 + 1;
#line 443 "csytrf_aa.f"
		    i__4 = jb + 1;
#line 443 "csytrf_aa.f"
		    cgemm_("No transpose", "Transpose", &i__3, &nj, &i__4, &
			    c_b19, &work[j3 - j1 + 1 + k1 * *n], n, &a[j2 + (
			    j1 - k2) * a_dim1], lda, &c_b15, &a[j3 + j2 * 
			    a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 448 "csytrf_aa.f"
		}

/*              Recover T( J+1, J ) */

#line 452 "csytrf_aa.f"
		i__1 = j + 1 + j * a_dim1;
#line 452 "csytrf_aa.f"
		a[i__1].r = alpha.r, a[i__1].i = alpha.i;
#line 453 "csytrf_aa.f"
	    }

/*           WORK(J+1, 1) stores H(J+1, 1) */

#line 457 "csytrf_aa.f"
	    i__1 = *n - j;
#line 457 "csytrf_aa.f"
	    ccopy_(&i__1, &a[j + 1 + (j + 1) * a_dim1], &c__1, &work[1], &
		    c__1);
#line 458 "csytrf_aa.f"
	}
#line 459 "csytrf_aa.f"
	goto L11;
#line 460 "csytrf_aa.f"
    }

#line 462 "csytrf_aa.f"
L20:
#line 463 "csytrf_aa.f"
    return 0;

/*     End of CSYTRF_AA */

} /* csytrf_aa__ */

