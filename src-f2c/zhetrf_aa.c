#line 1 "zhetrf_aa.f"
/* zhetrf_aa.f -- translated by f2c (version 20100827).
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

#line 1 "zhetrf_aa.f"
/* Table of constant values */

static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZHETRF_AA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRF_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrf_
aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrf_
aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrf_
aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER    UPLO */
/*       INTEGER      N, LDA, LWORK, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER      IPIV( * ) */
/*       COMPLEX*16   A( LDA, * ), WORK( * ) */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRF_AA computes the factorization of a complex hermitian matrix A */
/* > using the Aasen's algorithm.  The form of the factorization is */
/* > */
/* >    A = U*T*U**H  or  A = L*T*L**H */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is a hermitian tridiagonal matrix. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the hermitian matrix A.  If UPLO = 'U', the leading */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complex16HEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhetrf_aa__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int zlahef_aa__(char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, ftnlen);
    static integer j1, k1, k2, j2, j3, jb, nb, mj, nj;
    static doublecomplex alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), xerbla_(char *, integer *,
	     ftnlen);
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

#line 181 "zhetrf_aa.f"
    /* Parameter adjustments */
#line 181 "zhetrf_aa.f"
    a_dim1 = *lda;
#line 181 "zhetrf_aa.f"
    a_offset = 1 + a_dim1;
#line 181 "zhetrf_aa.f"
    a -= a_offset;
#line 181 "zhetrf_aa.f"
    --ipiv;
#line 181 "zhetrf_aa.f"
    --work;
#line 181 "zhetrf_aa.f"

#line 181 "zhetrf_aa.f"
    /* Function Body */
#line 181 "zhetrf_aa.f"
    nb = ilaenv_(&c__1, "ZHETRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

/*     Test the input parameters. */

#line 185 "zhetrf_aa.f"
    *info = 0;
#line 186 "zhetrf_aa.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 187 "zhetrf_aa.f"
    lquery = *lwork == -1;
#line 188 "zhetrf_aa.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 189 "zhetrf_aa.f"
	*info = -1;
#line 190 "zhetrf_aa.f"
    } else if (*n < 0) {
#line 191 "zhetrf_aa.f"
	*info = -2;
#line 192 "zhetrf_aa.f"
    } else if (*lda < max(1,*n)) {
#line 193 "zhetrf_aa.f"
	*info = -4;
#line 194 "zhetrf_aa.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 194 "zhetrf_aa.f"
	i__1 = 1, i__2 = *n << 1;
#line 194 "zhetrf_aa.f"
	if (*lwork < max(i__1,i__2) && ! lquery) {
#line 195 "zhetrf_aa.f"
	    *info = -7;
#line 196 "zhetrf_aa.f"
	}
#line 196 "zhetrf_aa.f"
    }

#line 198 "zhetrf_aa.f"
    if (*info == 0) {
#line 199 "zhetrf_aa.f"
	lwkopt = (nb + 1) * *n;
#line 200 "zhetrf_aa.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 201 "zhetrf_aa.f"
    }

#line 203 "zhetrf_aa.f"
    if (*info != 0) {
#line 204 "zhetrf_aa.f"
	i__1 = -(*info);
#line 204 "zhetrf_aa.f"
	xerbla_("ZHETRF_AA", &i__1, (ftnlen)9);
#line 205 "zhetrf_aa.f"
	return 0;
#line 206 "zhetrf_aa.f"
    } else if (lquery) {
#line 207 "zhetrf_aa.f"
	return 0;
#line 208 "zhetrf_aa.f"
    }

/*     Quick return */

#line 212 "zhetrf_aa.f"
    if (*n == 0) {
#line 213 "zhetrf_aa.f"
	return 0;
#line 214 "zhetrf_aa.f"
    }
#line 215 "zhetrf_aa.f"
    ipiv[1] = 1;
#line 216 "zhetrf_aa.f"
    if (*n == 1) {
#line 217 "zhetrf_aa.f"
	i__1 = a_dim1 + 1;
#line 217 "zhetrf_aa.f"
	i__2 = a_dim1 + 1;
#line 217 "zhetrf_aa.f"
	d__1 = a[i__2].r;
#line 217 "zhetrf_aa.f"
	a[i__1].r = d__1, a[i__1].i = 0.;
#line 218 "zhetrf_aa.f"
	i__1 = a_dim1 + 1;
#line 218 "zhetrf_aa.f"
	if (a[i__1].r == 0. && a[i__1].i == 0.) {
#line 219 "zhetrf_aa.f"
	    *info = 1;
#line 220 "zhetrf_aa.f"
	}
#line 221 "zhetrf_aa.f"
	return 0;
#line 222 "zhetrf_aa.f"
    }

/*     Adjubst block size based on the workspace size */

#line 226 "zhetrf_aa.f"
    if (*lwork < (nb + 1) * *n) {
#line 227 "zhetrf_aa.f"
	nb = (*lwork - *n) / *n;
#line 228 "zhetrf_aa.f"
    }

#line 230 "zhetrf_aa.f"
    if (upper) {

/*        ..................................................... */
/*        Factorize A as L*D*L**H using the upper triangle of A */
/*        ..................................................... */

/*        copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N)) */

#line 238 "zhetrf_aa.f"
	zcopy_(n, &a[a_dim1 + 1], lda, &work[1], &c__1);

/*        J is the main loop index, increasing from 1 to N in steps of */
/*        JB, where JB is the number of columns factorized by ZLAHEF; */
/*        JB is either NB, or N-J+1 for the last block */

#line 244 "zhetrf_aa.f"
	j = 0;
#line 245 "zhetrf_aa.f"
L10:
#line 246 "zhetrf_aa.f"
	if (j >= *n) {
#line 246 "zhetrf_aa.f"
	    goto L20;
#line 246 "zhetrf_aa.f"
	}

/*        each step of the main loop */
/*         J is the last column of the previous panel */
/*         J1 is the first column of the current panel */
/*         K1 identifies if the previous column of the panel has been */
/*          explicitly stored, e.g., K1=1 for the first panel, and */
/*          K1=0 for the rest */

#line 256 "zhetrf_aa.f"
	j1 = j + 1;
/* Computing MIN */
#line 257 "zhetrf_aa.f"
	i__1 = *n - j1 + 1;
#line 257 "zhetrf_aa.f"
	jb = min(i__1,nb);
#line 258 "zhetrf_aa.f"
	k1 = max(1,j) - j;

/*        Panel factorization */

#line 262 "zhetrf_aa.f"
	i__1 = 2 - k1;
#line 262 "zhetrf_aa.f"
	i__2 = *n - j;
#line 262 "zhetrf_aa.f"
	zlahef_aa__(uplo, &i__1, &i__2, &jb, &a[max(1,j) + (j + 1) * a_dim1], 
		lda, &ipiv[j + 1], &work[1], n, &work[*n * nb + 1], &iinfo, (
		ftnlen)1);
#line 266 "zhetrf_aa.f"
	if (iinfo > 0 && *info == 0) {
#line 267 "zhetrf_aa.f"
	    *info = iinfo + j;
#line 268 "zhetrf_aa.f"
	}

/*        Ajust IPIV and apply it back (J-th step picks (J+1)-th pivot) */

/* Computing MIN */
#line 272 "zhetrf_aa.f"
	i__2 = *n, i__3 = j + jb + 1;
#line 272 "zhetrf_aa.f"
	i__1 = min(i__2,i__3);
#line 272 "zhetrf_aa.f"
	for (j2 = j + 2; j2 <= i__1; ++j2) {
#line 273 "zhetrf_aa.f"
	    ipiv[j2] += j;
#line 274 "zhetrf_aa.f"
	    if (j2 != ipiv[j2] && j1 - k1 > 2) {
#line 275 "zhetrf_aa.f"
		i__2 = j1 - k1 - 2;
#line 275 "zhetrf_aa.f"
		zswap_(&i__2, &a[j2 * a_dim1 + 1], &c__1, &a[ipiv[j2] * 
			a_dim1 + 1], &c__1);
#line 277 "zhetrf_aa.f"
	    }
#line 278 "zhetrf_aa.f"
	}
#line 279 "zhetrf_aa.f"
	j += jb;

/*        Trailing submatrix update, where */
/*         the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and */
/*         WORK stores the current block of the auxiriarly matrix H */

#line 285 "zhetrf_aa.f"
	if (j < *n) {

/*          if the first panel and JB=1 (NB=1), then nothing to do */

#line 289 "zhetrf_aa.f"
	    if (j1 > 1 || jb > 1) {

/*              Merge rank-1 update with BLAS-3 update */

#line 293 "zhetrf_aa.f"
		d_cnjg(&z__1, &a[j + (j + 1) * a_dim1]);
#line 293 "zhetrf_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 294 "zhetrf_aa.f"
		i__1 = j + (j + 1) * a_dim1;
#line 294 "zhetrf_aa.f"
		a[i__1].r = 1., a[i__1].i = 0.;
#line 295 "zhetrf_aa.f"
		i__1 = *n - j;
#line 295 "zhetrf_aa.f"
		zcopy_(&i__1, &a[j - 1 + (j + 1) * a_dim1], lda, &work[j + 1 
			- j1 + 1 + jb * *n], &c__1);
#line 297 "zhetrf_aa.f"
		i__1 = *n - j;
#line 297 "zhetrf_aa.f"
		zscal_(&i__1, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);

/*              K1 identifies if the previous column of the panel has been */
/*               explicitly stored, e.g., K1=0 and K2=1 for the first panel, */
/*               and K1=1 and K2=0 for the rest */

#line 303 "zhetrf_aa.f"
		if (j1 > 1) {

/*                 Not first panel */

#line 307 "zhetrf_aa.f"
		    k2 = 1;
#line 308 "zhetrf_aa.f"
		} else {

/*                 First panel */

#line 312 "zhetrf_aa.f"
		    k2 = 0;

/*                 First update skips the first column */

#line 316 "zhetrf_aa.f"
		    --jb;
#line 317 "zhetrf_aa.f"
		}

#line 319 "zhetrf_aa.f"
		i__1 = *n;
#line 319 "zhetrf_aa.f"
		i__2 = nb;
#line 319 "zhetrf_aa.f"
		for (j2 = j + 1; i__2 < 0 ? j2 >= i__1 : j2 <= i__1; j2 += 
			i__2) {
/* Computing MIN */
#line 320 "zhetrf_aa.f"
		    i__3 = nb, i__4 = *n - j2 + 1;
#line 320 "zhetrf_aa.f"
		    nj = min(i__3,i__4);

/*                 Update (J2, J2) diagonal block with ZGEMV */

#line 324 "zhetrf_aa.f"
		    j3 = j2;
#line 325 "zhetrf_aa.f"
		    for (mj = nj - 1; mj >= 1; --mj) {
#line 326 "zhetrf_aa.f"
			i__3 = jb + 1;
#line 326 "zhetrf_aa.f"
			z__1.r = -1., z__1.i = -0.;
#line 326 "zhetrf_aa.f"
			zgemm_("Conjugate transpose", "Transpose", &c__1, &mj,
				 &i__3, &z__1, &a[j1 - k2 + j3 * a_dim1], lda,
				 &work[j3 - j1 + 1 + k1 * *n], n, &c_b2, &a[
				j3 + j3 * a_dim1], lda, (ftnlen)19, (ftnlen)9)
				;
#line 331 "zhetrf_aa.f"
			++j3;
#line 332 "zhetrf_aa.f"
		    }

/*                 Update off-diagonal block of J2-th block row with ZGEMM */

#line 336 "zhetrf_aa.f"
		    i__3 = *n - j3 + 1;
#line 336 "zhetrf_aa.f"
		    i__4 = jb + 1;
#line 336 "zhetrf_aa.f"
		    z__1.r = -1., z__1.i = -0.;
#line 336 "zhetrf_aa.f"
		    zgemm_("Conjugate transpose", "Transpose", &nj, &i__3, &
			    i__4, &z__1, &a[j1 - k2 + j2 * a_dim1], lda, &
			    work[j3 - j1 + 1 + k1 * *n], n, &c_b2, &a[j2 + j3 
			    * a_dim1], lda, (ftnlen)19, (ftnlen)9);
#line 341 "zhetrf_aa.f"
		}

/*              Recover T( J, J+1 ) */

#line 345 "zhetrf_aa.f"
		i__2 = j + (j + 1) * a_dim1;
#line 345 "zhetrf_aa.f"
		d_cnjg(&z__1, &alpha);
#line 345 "zhetrf_aa.f"
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 346 "zhetrf_aa.f"
	    }

/*           WORK(J+1, 1) stores H(J+1, 1) */

#line 350 "zhetrf_aa.f"
	    i__2 = *n - j;
#line 350 "zhetrf_aa.f"
	    zcopy_(&i__2, &a[j + 1 + (j + 1) * a_dim1], lda, &work[1], &c__1);
#line 351 "zhetrf_aa.f"
	}
#line 352 "zhetrf_aa.f"
	goto L10;
#line 353 "zhetrf_aa.f"
    } else {

/*        ..................................................... */
/*        Factorize A as L*D*L**H using the lower triangle of A */
/*        ..................................................... */

/*        copy first column A(1:N, 1) into H(1:N, 1) */
/*         (stored in WORK(1:N)) */

#line 362 "zhetrf_aa.f"
	zcopy_(n, &a[a_dim1 + 1], &c__1, &work[1], &c__1);

/*        J is the main loop index, increasing from 1 to N in steps of */
/*        JB, where JB is the number of columns factorized by ZLAHEF; */
/*        JB is either NB, or N-J+1 for the last block */

#line 368 "zhetrf_aa.f"
	j = 0;
#line 369 "zhetrf_aa.f"
L11:
#line 370 "zhetrf_aa.f"
	if (j >= *n) {
#line 370 "zhetrf_aa.f"
	    goto L20;
#line 370 "zhetrf_aa.f"
	}

/*        each step of the main loop */
/*         J is the last column of the previous panel */
/*         J1 is the first column of the current panel */
/*         K1 identifies if the previous column of the panel has been */
/*          explicitly stored, e.g., K1=1 for the first panel, and */
/*          K1=0 for the rest */

#line 380 "zhetrf_aa.f"
	j1 = j + 1;
/* Computing MIN */
#line 381 "zhetrf_aa.f"
	i__2 = *n - j1 + 1;
#line 381 "zhetrf_aa.f"
	jb = min(i__2,nb);
#line 382 "zhetrf_aa.f"
	k1 = max(1,j) - j;

/*        Panel factorization */

#line 386 "zhetrf_aa.f"
	i__2 = 2 - k1;
#line 386 "zhetrf_aa.f"
	i__1 = *n - j;
#line 386 "zhetrf_aa.f"
	zlahef_aa__(uplo, &i__2, &i__1, &jb, &a[j + 1 + max(1,j) * a_dim1], 
		lda, &ipiv[j + 1], &work[1], n, &work[*n * nb + 1], &iinfo, (
		ftnlen)1);
#line 389 "zhetrf_aa.f"
	if (iinfo > 0 && *info == 0) {
#line 390 "zhetrf_aa.f"
	    *info = iinfo + j;
#line 391 "zhetrf_aa.f"
	}

/*        Ajust IPIV and apply it back (J-th step picks (J+1)-th pivot) */

/* Computing MIN */
#line 395 "zhetrf_aa.f"
	i__1 = *n, i__3 = j + jb + 1;
#line 395 "zhetrf_aa.f"
	i__2 = min(i__1,i__3);
#line 395 "zhetrf_aa.f"
	for (j2 = j + 2; j2 <= i__2; ++j2) {
#line 396 "zhetrf_aa.f"
	    ipiv[j2] += j;
#line 397 "zhetrf_aa.f"
	    if (j2 != ipiv[j2] && j1 - k1 > 2) {
#line 398 "zhetrf_aa.f"
		i__1 = j1 - k1 - 2;
#line 398 "zhetrf_aa.f"
		zswap_(&i__1, &a[j2 + a_dim1], lda, &a[ipiv[j2] + a_dim1], 
			lda);
#line 400 "zhetrf_aa.f"
	    }
#line 401 "zhetrf_aa.f"
	}
#line 402 "zhetrf_aa.f"
	j += jb;

/*        Trailing submatrix update, where */
/*          A(J2+1, J1-1) stores L(J2+1, J1) and */
/*          WORK(J2+1, 1) stores H(J2+1, 1) */

#line 408 "zhetrf_aa.f"
	if (j < *n) {

/*          if the first panel and JB=1 (NB=1), then nothing to do */

#line 412 "zhetrf_aa.f"
	    if (j1 > 1 || jb > 1) {

/*              Merge rank-1 update with BLAS-3 update */

#line 416 "zhetrf_aa.f"
		d_cnjg(&z__1, &a[j + 1 + j * a_dim1]);
#line 416 "zhetrf_aa.f"
		alpha.r = z__1.r, alpha.i = z__1.i;
#line 417 "zhetrf_aa.f"
		i__2 = j + 1 + j * a_dim1;
#line 417 "zhetrf_aa.f"
		a[i__2].r = 1., a[i__2].i = 0.;
#line 418 "zhetrf_aa.f"
		i__2 = *n - j;
#line 418 "zhetrf_aa.f"
		zcopy_(&i__2, &a[j + 1 + (j - 1) * a_dim1], &c__1, &work[j + 
			1 - j1 + 1 + jb * *n], &c__1);
#line 420 "zhetrf_aa.f"
		i__2 = *n - j;
#line 420 "zhetrf_aa.f"
		zscal_(&i__2, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);

/*              K1 identifies if the previous column of the panel has been */
/*               explicitly stored, e.g., K1=0 and K2=1 for the first panel, */
/*               and K1=1 and K2=0 for the rest */

#line 426 "zhetrf_aa.f"
		if (j1 > 1) {

/*                 Not first panel */

#line 430 "zhetrf_aa.f"
		    k2 = 1;
#line 431 "zhetrf_aa.f"
		} else {

/*                 First panel */

#line 435 "zhetrf_aa.f"
		    k2 = 0;

/*                 First update skips the first column */

#line 439 "zhetrf_aa.f"
		    --jb;
#line 440 "zhetrf_aa.f"
		}

#line 442 "zhetrf_aa.f"
		i__2 = *n;
#line 442 "zhetrf_aa.f"
		i__1 = nb;
#line 442 "zhetrf_aa.f"
		for (j2 = j + 1; i__1 < 0 ? j2 >= i__2 : j2 <= i__2; j2 += 
			i__1) {
/* Computing MIN */
#line 443 "zhetrf_aa.f"
		    i__3 = nb, i__4 = *n - j2 + 1;
#line 443 "zhetrf_aa.f"
		    nj = min(i__3,i__4);

/*                 Update (J2, J2) diagonal block with ZGEMV */

#line 447 "zhetrf_aa.f"
		    j3 = j2;
#line 448 "zhetrf_aa.f"
		    for (mj = nj - 1; mj >= 1; --mj) {
#line 449 "zhetrf_aa.f"
			i__3 = jb + 1;
#line 449 "zhetrf_aa.f"
			z__1.r = -1., z__1.i = -0.;
#line 449 "zhetrf_aa.f"
			zgemm_("No transpose", "Conjugate transpose", &mj, &
				c__1, &i__3, &z__1, &work[j3 - j1 + 1 + k1 * *
				n], n, &a[j3 + (j1 - k2) * a_dim1], lda, &
				c_b2, &a[j3 + j3 * a_dim1], lda, (ftnlen)12, (
				ftnlen)19);
#line 454 "zhetrf_aa.f"
			++j3;
#line 455 "zhetrf_aa.f"
		    }

/*                 Update off-diagonal block of J2-th block column with ZGEMM */

#line 459 "zhetrf_aa.f"
		    i__3 = *n - j3 + 1;
#line 459 "zhetrf_aa.f"
		    i__4 = jb + 1;
#line 459 "zhetrf_aa.f"
		    z__1.r = -1., z__1.i = -0.;
#line 459 "zhetrf_aa.f"
		    zgemm_("No transpose", "Conjugate transpose", &i__3, &nj, 
			    &i__4, &z__1, &work[j3 - j1 + 1 + k1 * *n], n, &a[
			    j2 + (j1 - k2) * a_dim1], lda, &c_b2, &a[j3 + j2 *
			     a_dim1], lda, (ftnlen)12, (ftnlen)19);
#line 464 "zhetrf_aa.f"
		}

/*              Recover T( J+1, J ) */

#line 468 "zhetrf_aa.f"
		i__1 = j + 1 + j * a_dim1;
#line 468 "zhetrf_aa.f"
		d_cnjg(&z__1, &alpha);
#line 468 "zhetrf_aa.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 469 "zhetrf_aa.f"
	    }

/*           WORK(J+1, 1) stores H(J+1, 1) */

#line 473 "zhetrf_aa.f"
	    i__1 = *n - j;
#line 473 "zhetrf_aa.f"
	    zcopy_(&i__1, &a[j + 1 + (j + 1) * a_dim1], &c__1, &work[1], &
		    c__1);
#line 474 "zhetrf_aa.f"
	}
#line 475 "zhetrf_aa.f"
	goto L11;
#line 476 "zhetrf_aa.f"
    }

#line 478 "zhetrf_aa.f"
L20:
#line 479 "zhetrf_aa.f"
    return 0;

/*     End of ZHETRF_AA */

} /* zhetrf_aa__ */

