#line 1 "zhetrs_3.f"
/* zhetrs_3.f -- translated by f2c (version 20100827).
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

#line 1 "zhetrs_3.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b ZHETRS_3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRS_3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_
3.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_
3.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_
3.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, */
/*                            INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > ZHETRS_3 solves a system of linear equations A * X = B with a complex */
/* > Hermitian matrix A using the factorization computed */
/* > by ZHETRF_RK or ZHETRF_BK: */
/* > */
/* >    A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is Hermitian and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This algorithm is using Level 3 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix: */
/* >          = 'U':  Upper triangular, form is A = P*U*D*(U**H)*(P**T); */
/* >          = 'L':  Lower triangular, form is A = P*L*D*(L**H)*(P**T). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          Diagonal of the block diagonal matrix D and factors U or L */
/* >          as computed by ZHETRF_RK and ZHETRF_BK: */
/* >            a) ONLY diagonal elements of the Hermitian block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                should be provided on entry in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is COMPLEX*16 array, dimension (N) */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the Hermitian block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          NOTE: For 1-by-1 diagonal block D(k), where */
/* >          1 <= k <= N, the element E(k) is not referenced in both */
/* >          UPLO = 'U' or UPLO = 'L' cases. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZHETRF_RK or ZHETRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side matrix B. */
/* >          On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup complex16HEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  June 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zhetrs_3__(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	doublecomplex *b, integer *ldb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex ak, bk;
    static integer kp;
    static doublecomplex akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex denom;
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ztrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), zdscal_(integer *, doublereal 
	    *, doublecomplex *, integer *);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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
/*     .. Executable Statements .. */

#line 206 "zhetrs_3.f"
    /* Parameter adjustments */
#line 206 "zhetrs_3.f"
    a_dim1 = *lda;
#line 206 "zhetrs_3.f"
    a_offset = 1 + a_dim1;
#line 206 "zhetrs_3.f"
    a -= a_offset;
#line 206 "zhetrs_3.f"
    --e;
#line 206 "zhetrs_3.f"
    --ipiv;
#line 206 "zhetrs_3.f"
    b_dim1 = *ldb;
#line 206 "zhetrs_3.f"
    b_offset = 1 + b_dim1;
#line 206 "zhetrs_3.f"
    b -= b_offset;
#line 206 "zhetrs_3.f"

#line 206 "zhetrs_3.f"
    /* Function Body */
#line 206 "zhetrs_3.f"
    *info = 0;
#line 207 "zhetrs_3.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 208 "zhetrs_3.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 209 "zhetrs_3.f"
	*info = -1;
#line 210 "zhetrs_3.f"
    } else if (*n < 0) {
#line 211 "zhetrs_3.f"
	*info = -2;
#line 212 "zhetrs_3.f"
    } else if (*nrhs < 0) {
#line 213 "zhetrs_3.f"
	*info = -3;
#line 214 "zhetrs_3.f"
    } else if (*lda < max(1,*n)) {
#line 215 "zhetrs_3.f"
	*info = -5;
#line 216 "zhetrs_3.f"
    } else if (*ldb < max(1,*n)) {
#line 217 "zhetrs_3.f"
	*info = -9;
#line 218 "zhetrs_3.f"
    }
#line 219 "zhetrs_3.f"
    if (*info != 0) {
#line 220 "zhetrs_3.f"
	i__1 = -(*info);
#line 220 "zhetrs_3.f"
	xerbla_("ZHETRS_3", &i__1, (ftnlen)8);
#line 221 "zhetrs_3.f"
	return 0;
#line 222 "zhetrs_3.f"
    }

/*     Quick return if possible */

#line 226 "zhetrs_3.f"
    if (*n == 0 || *nrhs == 0) {
#line 226 "zhetrs_3.f"
	return 0;
#line 226 "zhetrs_3.f"
    }

#line 229 "zhetrs_3.f"
    if (upper) {

/*        Begin Upper */

/*        Solve A*X = B, where A = U*D*U**H. */

/*        P**T * B */

/*        Interchange rows K and IPIV(K) of matrix B in the same order */
/*        that the formation order of IPIV(I) vector for Upper case. */

/*        (We can do the simple loop over IPIV with decrement -1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 244 "zhetrs_3.f"
	for (k = *n; k >= 1; --k) {
#line 245 "zhetrs_3.f"
	    kp = (i__1 = ipiv[k], abs(i__1));
#line 246 "zhetrs_3.f"
	    if (kp != k) {
#line 247 "zhetrs_3.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 248 "zhetrs_3.f"
	    }
#line 249 "zhetrs_3.f"
	}

/*        Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 253 "zhetrs_3.f"
	ztrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 257 "zhetrs_3.f"
	i__ = *n;
#line 258 "zhetrs_3.f"
	while(i__ >= 1) {
#line 259 "zhetrs_3.f"
	    if (ipiv[i__] > 0) {
#line 260 "zhetrs_3.f"
		i__1 = i__ + i__ * a_dim1;
#line 260 "zhetrs_3.f"
		s = 1. / a[i__1].r;
#line 261 "zhetrs_3.f"
		zdscal_(nrhs, &s, &b[i__ + b_dim1], ldb);
#line 262 "zhetrs_3.f"
	    } else if (i__ > 1) {
#line 263 "zhetrs_3.f"
		i__1 = i__;
#line 263 "zhetrs_3.f"
		akm1k.r = e[i__1].r, akm1k.i = e[i__1].i;
#line 264 "zhetrs_3.f"
		z_div(&z__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
#line 264 "zhetrs_3.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 265 "zhetrs_3.f"
		d_cnjg(&z__2, &akm1k);
#line 265 "zhetrs_3.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &z__2);
#line 265 "zhetrs_3.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 266 "zhetrs_3.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 266 "zhetrs_3.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 266 "zhetrs_3.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 267 "zhetrs_3.f"
		i__1 = *nrhs;
#line 267 "zhetrs_3.f"
		for (j = 1; j <= i__1; ++j) {
#line 268 "zhetrs_3.f"
		    z_div(&z__1, &b[i__ - 1 + j * b_dim1], &akm1k);
#line 268 "zhetrs_3.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 269 "zhetrs_3.f"
		    d_cnjg(&z__2, &akm1k);
#line 269 "zhetrs_3.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &z__2);
#line 269 "zhetrs_3.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 270 "zhetrs_3.f"
		    i__2 = i__ - 1 + j * b_dim1;
#line 270 "zhetrs_3.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 270 "zhetrs_3.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 270 "zhetrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 270 "zhetrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 271 "zhetrs_3.f"
		    i__2 = i__ + j * b_dim1;
#line 271 "zhetrs_3.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 271 "zhetrs_3.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 271 "zhetrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 271 "zhetrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 272 "zhetrs_3.f"
		}
#line 273 "zhetrs_3.f"
		--i__;
#line 274 "zhetrs_3.f"
	    }
#line 275 "zhetrs_3.f"
	    --i__;
#line 276 "zhetrs_3.f"
	}

/*        Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ] */

#line 280 "zhetrs_3.f"
	ztrsm_("L", "U", "C", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ] */

/*        Interchange rows K and IPIV(K) of matrix B in reverse order */
/*        from the formation order of IPIV(I) vector for Upper case. */

/*        (We can do the simple loop over IPIV with increment 1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 291 "zhetrs_3.f"
	i__1 = *n;
#line 291 "zhetrs_3.f"
	for (k = 1; k <= i__1; ++k) {
#line 292 "zhetrs_3.f"
	    kp = (i__2 = ipiv[k], abs(i__2));
#line 293 "zhetrs_3.f"
	    if (kp != k) {
#line 294 "zhetrs_3.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 295 "zhetrs_3.f"
	    }
#line 296 "zhetrs_3.f"
	}

#line 298 "zhetrs_3.f"
    } else {

/*        Begin Lower */

/*        Solve A*X = B, where A = L*D*L**H. */

/*        P**T * B */
/*        Interchange rows K and IPIV(K) of matrix B in the same order */
/*        that the formation order of IPIV(I) vector for Lower case. */

/*        (We can do the simple loop over IPIV with increment 1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 312 "zhetrs_3.f"
	i__1 = *n;
#line 312 "zhetrs_3.f"
	for (k = 1; k <= i__1; ++k) {
#line 313 "zhetrs_3.f"
	    kp = (i__2 = ipiv[k], abs(i__2));
#line 314 "zhetrs_3.f"
	    if (kp != k) {
#line 315 "zhetrs_3.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 316 "zhetrs_3.f"
	    }
#line 317 "zhetrs_3.f"
	}

/*        Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 321 "zhetrs_3.f"
	ztrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 325 "zhetrs_3.f"
	i__ = 1;
#line 326 "zhetrs_3.f"
	while(i__ <= *n) {
#line 327 "zhetrs_3.f"
	    if (ipiv[i__] > 0) {
#line 328 "zhetrs_3.f"
		i__1 = i__ + i__ * a_dim1;
#line 328 "zhetrs_3.f"
		s = 1. / a[i__1].r;
#line 329 "zhetrs_3.f"
		zdscal_(nrhs, &s, &b[i__ + b_dim1], ldb);
#line 330 "zhetrs_3.f"
	    } else if (i__ < *n) {
#line 331 "zhetrs_3.f"
		i__1 = i__;
#line 331 "zhetrs_3.f"
		akm1k.r = e[i__1].r, akm1k.i = e[i__1].i;
#line 332 "zhetrs_3.f"
		d_cnjg(&z__2, &akm1k);
#line 332 "zhetrs_3.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &z__2);
#line 332 "zhetrs_3.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 333 "zhetrs_3.f"
		z_div(&z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
#line 333 "zhetrs_3.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 334 "zhetrs_3.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 334 "zhetrs_3.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 334 "zhetrs_3.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 335 "zhetrs_3.f"
		i__1 = *nrhs;
#line 335 "zhetrs_3.f"
		for (j = 1; j <= i__1; ++j) {
#line 336 "zhetrs_3.f"
		    d_cnjg(&z__2, &akm1k);
#line 336 "zhetrs_3.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &z__2);
#line 336 "zhetrs_3.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 337 "zhetrs_3.f"
		    z_div(&z__1, &b[i__ + 1 + j * b_dim1], &akm1k);
#line 337 "zhetrs_3.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 338 "zhetrs_3.f"
		    i__2 = i__ + j * b_dim1;
#line 338 "zhetrs_3.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 338 "zhetrs_3.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 338 "zhetrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 338 "zhetrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 339 "zhetrs_3.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 339 "zhetrs_3.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 339 "zhetrs_3.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 339 "zhetrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 339 "zhetrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 340 "zhetrs_3.f"
		}
#line 341 "zhetrs_3.f"
		++i__;
#line 342 "zhetrs_3.f"
	    }
#line 343 "zhetrs_3.f"
	    ++i__;
#line 344 "zhetrs_3.f"
	}

/*        Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ] */

#line 348 "zhetrs_3.f"
	ztrsm_("L", "L", "C", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ] */

/*        Interchange rows K and IPIV(K) of matrix B in reverse order */
/*        from the formation order of IPIV(I) vector for Lower case. */

/*        (We can do the simple loop over IPIV with decrement -1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 359 "zhetrs_3.f"
	for (k = *n; k >= 1; --k) {
#line 360 "zhetrs_3.f"
	    kp = (i__1 = ipiv[k], abs(i__1));
#line 361 "zhetrs_3.f"
	    if (kp != k) {
#line 362 "zhetrs_3.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 363 "zhetrs_3.f"
	    }
#line 364 "zhetrs_3.f"
	}

/*        END Lower */

#line 368 "zhetrs_3.f"
    }

#line 370 "zhetrs_3.f"
    return 0;

/*     End of ZHETRS_3 */

} /* zhetrs_3__ */

