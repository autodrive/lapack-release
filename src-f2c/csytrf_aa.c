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
static doublecomplex c_b16 = {1.,0.};
static doublecomplex c_b20 = {-1.,-0.};

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
	     integer *, doublecomplex *, integer *, ftnlen);
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
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
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

#line 181 "csytrf_aa.f"
    /* Parameter adjustments */
#line 181 "csytrf_aa.f"
    a_dim1 = *lda;
#line 181 "csytrf_aa.f"
    a_offset = 1 + a_dim1;
#line 181 "csytrf_aa.f"
    a -= a_offset;
#line 181 "csytrf_aa.f"
    --ipiv;
#line 181 "csytrf_aa.f"
    --work;
#line 181 "csytrf_aa.f"

#line 181 "csytrf_aa.f"
    /* Function Body */
#line 181 "csytrf_aa.f"
    nb = ilaenv_(&c__1, "CSYTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

/*     Test the input parameters. */

#line 185 "csytrf_aa.f"
    *info = 0;
#line 186 "csytrf_aa.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 187 "csytrf_aa.f"
    lquery = *lwork == -1;
#line 188 "csytrf_aa.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 189 "csytrf_aa.f"
	*info = -1;
#line 190 "csytrf_aa.f"
    } else if (*n < 0) {
#line 191 "csytrf_aa.f"
	*info = -2;
#line 192 "csytrf_aa.f"
    } else if (*lda < max(1,*n)) {
#line 193 "csytrf_aa.f"
	*info = -4;
#line 194 "csytrf_aa.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 194 "csytrf_aa.f"
	i__1 = 1, i__2 = *n << 1;
#line 194 "csytrf_aa.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 195 "csytrf_aa.f"
	    *info = -7;
#line 196 "csytrf_aa.f"
	}
#line 196 "csytrf_aa.f"
    }

#line 198 "csytrf_aa.f"
    if (*info == 0) {
#line 199 "csytrf_aa.f"
	lwkopt = (nb + 1) * *n;
#line 200 "csytrf_aa.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 201 "csytrf_aa.f"
    }

#line 203 "csytrf_aa.f"
    if (*info != 0) {
#line 204 "csytrf_aa.f"
	i__1 = -(*info);
#line 204 "csytrf_aa.f"
	xerbla_("CSYTRF_AA", &i__1, (ftnlen)9);
#line 205 "csytrf_aa.f"
	return 0;
#line 206 "csytrf_aa.f"
    } else if (lquery) {
#line 207 "csytrf_aa.f"
	return 0;
#line 208 "csytrf_aa.f"
    }

/*     Quick return */

#line 212 "csytrf_aa.f"
    if (*n == 0) {
#line 213 "csytrf_aa.f"
	return 0;
#line 214 "csytrf_aa.f"
    }
#line 215 "csytrf_aa.f"
    ipiv[1] = 1;
#line 216 "csytrf_aa.f"
    if (*n == 1) {
#line 217 "csytrf_aa.f"
	i__1 = a_dim1 + 1;
#line 217 "csytrf_aa.f"
	if (a[i__1].r == 0. && a[i__1].i == 0.) {
#line 218 "csytrf_aa.f"
	    *info = 1;
#line 219 "csytrf_aa.f"
	}
#line 220 "csytrf_aa.f"
	return 0;
#line 221 "csytrf_aa.f"
    }

/*     Adjubst block size based on the workspace size */

#line 225 "csytrf_aa.f"
    if (*lwork < (nb + 1) * *n) {
#line 226 "csytrf_aa.f"
	nb = (*lwork - *n) / *n;
#line 227 "csytrf_aa.f"
    }

#line 229 "csytrf_aa.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the upper triangle of A */
/*        ..................................................... */

/*        Copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N)) */

#line 237 "csytrf_aa.f"
	ccopy_(n, &a[a_dim1 + 1], lda, &work[1], &c__1);

/*        J is the main loop index, increasing from 1 to N in steps of */
/*        JB, where JB is the number of columns factorized by CLASYF; */
/*        JB is either NB, or N-J+1 for the last block */

#line 243 "csytrf_aa.f"
	j = 0;
#line 244 "csytrf_aa.f"
L10:
#line 245 "csytrf_aa.f"
	if (j >= *n) {
#line 245 "csytrf_aa.f"
	    goto L20;
#line 245 "csytrf_aa.f"
	}

/*        each step of the main loop */
/*         J is the last column of the previous panel */
/*         J1 is the first column of the current panel */
/*         K1 identifies if the previous column of the panel has been */
/*          explicitly stored, e.g., K1=1 for the first panel, and */
/*          K1=0 for the rest */

#line 255 "csytrf_aa.f"
	j1 = j + 1;
/* Computing MIN */
#line 256 "csytrf_aa.f"
	i__1 = *n - j1 + 1;
#line 256 "csytrf_aa.f"
	jb = min(i__1,nb);
#line 257 "csytrf_aa.f"
	k1 = max(1,j) - j;

/*        Panel factorization */

#line 261 "csytrf_aa.f"
	i__1 = 2 - k1;
#line 261 "csytrf_aa.f"
	i__2 = *n - j;
#line 261 "csytrf_aa.f"
	clasyf_aa__(uplo, &i__1, &i__2, &jb, &a[max(1,j) + (j + 1) * a_dim1], 
		lda, &ipiv[j + 1], &work[1], n, &work[*n * nb + 1], &iinfo, (
		ftnlen)1);
#line 265 "csytrf_aa.f"
	if (iinfo > 0 && *info == 0) {
#line 266 "csytrf_aa.f"
	    *info = iinfo + j;
#line 267 "csytrf_aa.f"
	}

/*        Ajust IPIV and apply it back (J-th step picks (J+1)-th pivot) */

/* Computing MIN */
#line 271 "csytrf_aa.f"
	i__2 = *n, i__3 = j + jb + 1;
#line 271 "csytrf_aa.f"
	i__1 = min(i__2,i__3);
#line 271 "csytrf_aa.f"
	for (j2 = j + 2; j2 <= i__1; ++j2) {
#line 272 "csytrf_aa.f"
	    ipiv[j2] += j;
#line 273 "csytrf_aa.f"
	    if (j2 != ipiv[j2] && j1 - k1 > 2) {
#line 274 "csytrf_aa.f"
		i__2 = j1 - k1 - 2;
#line 274 "csytrf_aa.f"
		cswap_(&i__2, &a[j2 * a_dim1 + 1], &c__1, &a[ipiv[j2] * 
			a_dim1 + 1], &c__1);
#line 276 "csytrf_aa.f"
	    }
#line 277 "csytrf_aa.f"
	}
#line 278 "csytrf_aa.f"
	j += jb;

/*        Trailing submatrix update, where */
/*         the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and */
/*         WORK stores the current block of the auxiriarly matrix H */

#line 284 "csytrf_aa.f"
	if (j < *n) {

/*           If first panel and JB=1 (NB=1), then nothing to do */

#line 288 "csytrf_aa.f"
	    if (j1 > 1 || jb > 1) {

/*              Merge rank-1 update with BLAS-3 update */

#line 292 "csytrf_aa.f"
		i__1 = j + (j + 1) * a_dim1;
#line 292 "csytrf_aa.f"
		alpha.r = a[i__1].r, alpha.i = a[i__1].i;
#line 293 "csytrf_aa.f"
		i__1 = j + (j + 1) * a_dim1;
#line 293 "csytrf_aa.f"
		a[i__1].r = 1., a[i__1].i = 0.;
#line 294 "csytrf_aa.f"
		i__1 = *n - j;
#line 294 "csytrf_aa.f"
		ccopy_(&i__1, &a[j - 1 + (j + 1) * a_dim1], lda, &work[j + 1 
			- j1 + 1 + jb * *n], &c__1);
#line 296 "csytrf_aa.f"
		i__1 = *n - j;
#line 296 "csytrf_aa.f"
		cscal_(&i__1, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);

/*              K1 identifies if the previous column of the panel has been */
/*               explicitly stored, e.g., K1=1 and K2= 0 for the first panel, */
/*               while K1=0 and K2=1 for the rest */

#line 302 "csytrf_aa.f"
		if (j1 > 1) {

/*                 Not first panel */

#line 306 "csytrf_aa.f"
		    k2 = 1;
#line 307 "csytrf_aa.f"
		} else {

/*                 First panel */

#line 311 "csytrf_aa.f"
		    k2 = 0;

/*                 First update skips the first column */

#line 315 "csytrf_aa.f"
		    --jb;
#line 316 "csytrf_aa.f"
		}

#line 318 "csytrf_aa.f"
		i__1 = *n;
#line 318 "csytrf_aa.f"
		i__2 = nb;
#line 318 "csytrf_aa.f"
		for (j2 = j + 1; i__2 < 0 ? j2 >= i__1 : j2 <= i__1; j2 += 
			i__2) {
/* Computing MIN */
#line 319 "csytrf_aa.f"
		    i__3 = nb, i__4 = *n - j2 + 1;
#line 319 "csytrf_aa.f"
		    nj = min(i__3,i__4);

/*                 Update (J2, J2) diagonal block with CGEMV */

#line 323 "csytrf_aa.f"
		    j3 = j2;
#line 324 "csytrf_aa.f"
		    for (mj = nj - 1; mj >= 1; --mj) {
#line 325 "csytrf_aa.f"
			i__3 = jb + 1;
#line 325 "csytrf_aa.f"
			cgemv_("No transpose", &mj, &i__3, &c_b20, &work[j3 - 
				j1 + 1 + k1 * *n], n, &a[j1 - k2 + j3 * 
				a_dim1], &c__1, &c_b16, &a[j3 + j3 * a_dim1], 
				lda, (ftnlen)12);
#line 329 "csytrf_aa.f"
			++j3;
#line 330 "csytrf_aa.f"
		    }

/*                 Update off-diagonal block of J2-th block row with CGEMM */

#line 334 "csytrf_aa.f"
		    i__3 = *n - j3 + 1;
#line 334 "csytrf_aa.f"
		    i__4 = jb + 1;
#line 334 "csytrf_aa.f"
		    cgemm_("Transpose", "Transpose", &nj, &i__3, &i__4, &
			    c_b20, &a[j1 - k2 + j2 * a_dim1], lda, &work[j3 - 
			    j1 + 1 + k1 * *n], n, &c_b16, &a[j2 + j3 * a_dim1]
			    , lda, (ftnlen)9, (ftnlen)9);
#line 339 "csytrf_aa.f"
		}

/*              Recover T( J, J+1 ) */

#line 343 "csytrf_aa.f"
		i__2 = j + (j + 1) * a_dim1;
#line 343 "csytrf_aa.f"
		a[i__2].r = alpha.r, a[i__2].i = alpha.i;
#line 344 "csytrf_aa.f"
	    }

/*           WORK(J+1, 1) stores H(J+1, 1) */

#line 348 "csytrf_aa.f"
	    i__2 = *n - j;
#line 348 "csytrf_aa.f"
	    ccopy_(&i__2, &a[j + 1 + (j + 1) * a_dim1], lda, &work[1], &c__1);
#line 349 "csytrf_aa.f"
	}
#line 350 "csytrf_aa.f"
	goto L10;
#line 351 "csytrf_aa.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**T using the lower triangle of A */
/*        ..................................................... */

/*        copy first column A(1:N, 1) into H(1:N, 1) */
/*         (stored in WORK(1:N)) */

#line 360 "csytrf_aa.f"
	ccopy_(n, &a[a_dim1 + 1], &c__1, &work[1], &c__1);

/*        J is the main loop index, increasing from 1 to N in steps of */
/*        JB, where JB is the number of columns factorized by CLASYF; */
/*        JB is either NB, or N-J+1 for the last block */

#line 366 "csytrf_aa.f"
	j = 0;
#line 367 "csytrf_aa.f"
L11:
#line 368 "csytrf_aa.f"
	if (j >= *n) {
#line 368 "csytrf_aa.f"
	    goto L20;
#line 368 "csytrf_aa.f"
	}

/*        each step of the main loop */
/*         J is the last column of the previous panel */
/*         J1 is the first column of the current panel */
/*         K1 identifies if the previous column of the panel has been */
/*          explicitly stored, e.g., K1=1 for the first panel, and */
/*          K1=0 for the rest */

#line 378 "csytrf_aa.f"
	j1 = j + 1;
/* Computing MIN */
#line 379 "csytrf_aa.f"
	i__2 = *n - j1 + 1;
#line 379 "csytrf_aa.f"
	jb = min(i__2,nb);
#line 380 "csytrf_aa.f"
	k1 = max(1,j) - j;

/*        Panel factorization */

#line 384 "csytrf_aa.f"
	i__2 = 2 - k1;
#line 384 "csytrf_aa.f"
	i__1 = *n - j;
#line 384 "csytrf_aa.f"
	clasyf_aa__(uplo, &i__2, &i__1, &jb, &a[j + 1 + max(1,j) * a_dim1], 
		lda, &ipiv[j + 1], &work[1], n, &work[*n * nb + 1], &iinfo, (
		ftnlen)1);
#line 387 "csytrf_aa.f"
	if (iinfo > 0 && *info == 0) {
#line 388 "csytrf_aa.f"
	    *info = iinfo + j;
#line 389 "csytrf_aa.f"
	}

/*        Ajust IPIV and apply it back (J-th step picks (J+1)-th pivot) */

/* Computing MIN */
#line 393 "csytrf_aa.f"
	i__1 = *n, i__3 = j + jb + 1;
#line 393 "csytrf_aa.f"
	i__2 = min(i__1,i__3);
#line 393 "csytrf_aa.f"
	for (j2 = j + 2; j2 <= i__2; ++j2) {
#line 394 "csytrf_aa.f"
	    ipiv[j2] += j;
#line 395 "csytrf_aa.f"
	    if (j2 != ipiv[j2] && j1 - k1 > 2) {
#line 396 "csytrf_aa.f"
		i__1 = j1 - k1 - 2;
#line 396 "csytrf_aa.f"
		cswap_(&i__1, &a[j2 + a_dim1], lda, &a[ipiv[j2] + a_dim1], 
			lda);
#line 398 "csytrf_aa.f"
	    }
#line 399 "csytrf_aa.f"
	}
#line 400 "csytrf_aa.f"
	j += jb;

/*        Trailing submatrix update, where */
/*          A(J2+1, J1-1) stores L(J2+1, J1) and */
/*          WORK(J2+1, 1) stores H(J2+1, 1) */

#line 406 "csytrf_aa.f"
	if (j < *n) {

/*           if first panel and JB=1 (NB=1), then nothing to do */

#line 410 "csytrf_aa.f"
	    if (j1 > 1 || jb > 1) {

/*              Merge rank-1 update with BLAS-3 update */

#line 414 "csytrf_aa.f"
		i__2 = j + 1 + j * a_dim1;
#line 414 "csytrf_aa.f"
		alpha.r = a[i__2].r, alpha.i = a[i__2].i;
#line 415 "csytrf_aa.f"
		i__2 = j + 1 + j * a_dim1;
#line 415 "csytrf_aa.f"
		a[i__2].r = 1., a[i__2].i = 0.;
#line 416 "csytrf_aa.f"
		i__2 = *n - j;
#line 416 "csytrf_aa.f"
		ccopy_(&i__2, &a[j + 1 + (j - 1) * a_dim1], &c__1, &work[j + 
			1 - j1 + 1 + jb * *n], &c__1);
#line 418 "csytrf_aa.f"
		i__2 = *n - j;
#line 418 "csytrf_aa.f"
		cscal_(&i__2, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);

/*              K1 identifies if the previous column of the panel has been */
/*               explicitly stored, e.g., K1=1 and K2= 0 for the first panel, */
/*               while K1=0 and K2=1 for the rest */

#line 424 "csytrf_aa.f"
		if (j1 > 1) {

/*                 Not first panel */

#line 428 "csytrf_aa.f"
		    k2 = 1;
#line 429 "csytrf_aa.f"
		} else {

/*                 First panel */

#line 433 "csytrf_aa.f"
		    k2 = 0;

/*                 First update skips the first column */

#line 437 "csytrf_aa.f"
		    --jb;
#line 438 "csytrf_aa.f"
		}

#line 440 "csytrf_aa.f"
		i__2 = *n;
#line 440 "csytrf_aa.f"
		i__1 = nb;
#line 440 "csytrf_aa.f"
		for (j2 = j + 1; i__1 < 0 ? j2 >= i__2 : j2 <= i__2; j2 += 
			i__1) {
/* Computing MIN */
#line 441 "csytrf_aa.f"
		    i__3 = nb, i__4 = *n - j2 + 1;
#line 441 "csytrf_aa.f"
		    nj = min(i__3,i__4);

/*                 Update (J2, J2) diagonal block with CGEMV */

#line 445 "csytrf_aa.f"
		    j3 = j2;
#line 446 "csytrf_aa.f"
		    for (mj = nj - 1; mj >= 1; --mj) {
#line 447 "csytrf_aa.f"
			i__3 = jb + 1;
#line 447 "csytrf_aa.f"
			cgemv_("No transpose", &mj, &i__3, &c_b20, &work[j3 - 
				j1 + 1 + k1 * *n], n, &a[j3 + (j1 - k2) * 
				a_dim1], lda, &c_b16, &a[j3 + j3 * a_dim1], &
				c__1, (ftnlen)12);
#line 451 "csytrf_aa.f"
			++j3;
#line 452 "csytrf_aa.f"
		    }

/*                 Update off-diagonal block in J2-th block column with CGEMM */

#line 456 "csytrf_aa.f"
		    i__3 = *n - j3 + 1;
#line 456 "csytrf_aa.f"
		    i__4 = jb + 1;
#line 456 "csytrf_aa.f"
		    cgemm_("No transpose", "Transpose", &i__3, &nj, &i__4, &
			    c_b20, &work[j3 - j1 + 1 + k1 * *n], n, &a[j2 + (
			    j1 - k2) * a_dim1], lda, &c_b16, &a[j3 + j2 * 
			    a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 461 "csytrf_aa.f"
		}

/*              Recover T( J+1, J ) */

#line 465 "csytrf_aa.f"
		i__1 = j + 1 + j * a_dim1;
#line 465 "csytrf_aa.f"
		a[i__1].r = alpha.r, a[i__1].i = alpha.i;
#line 466 "csytrf_aa.f"
	    }

/*           WORK(J+1, 1) stores H(J+1, 1) */

#line 470 "csytrf_aa.f"
	    i__1 = *n - j;
#line 470 "csytrf_aa.f"
	    ccopy_(&i__1, &a[j + 1 + (j + 1) * a_dim1], &c__1, &work[1], &
		    c__1);
#line 471 "csytrf_aa.f"
	}
#line 472 "csytrf_aa.f"
	goto L11;
#line 473 "csytrf_aa.f"
    }

#line 475 "csytrf_aa.f"
L20:
#line 476 "csytrf_aa.f"
    return 0;

/*     End of CSYTRF_AA */

} /* csytrf_aa__ */

