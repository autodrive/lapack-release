#line 1 "zsytrs_3.f"
/* zsytrs_3.f -- translated by f2c (version 20100827).
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

#line 1 "zsytrs_3.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b ZSYTRS_3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRS_3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs_
3.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs_
3.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs_
3.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, */
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
/* > ZSYTRS_3 solves a system of linear equations A * X = B with a complex */
/* > symmetric matrix A using the factorization computed */
/* > by ZSYTRF_RK or ZSYTRF_BK: */
/* > */
/* >    A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is symmetric and block */
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
/* >          = 'U':  Upper triangular, form is A = P*U*D*(U**T)*(P**T); */
/* >          = 'L':  Lower triangular, form is A = P*L*D*(L**T)*(P**T). */
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
/* >          as computed by ZSYTRF_RK and ZSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
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
/* >          elements of the symmetric block diagonal matrix D */
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
/* >          as determined by ZSYTRF_RK or ZSYTRF_BK. */
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

/* > \ingroup complex16SYcomputational */

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
/* Subroutine */ int zsytrs_3__(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	doublecomplex *b, integer *ldb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex ak, bk;
    static integer kp;
    static doublecomplex akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex denom;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ztrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

#line 205 "zsytrs_3.f"
    /* Parameter adjustments */
#line 205 "zsytrs_3.f"
    a_dim1 = *lda;
#line 205 "zsytrs_3.f"
    a_offset = 1 + a_dim1;
#line 205 "zsytrs_3.f"
    a -= a_offset;
#line 205 "zsytrs_3.f"
    --e;
#line 205 "zsytrs_3.f"
    --ipiv;
#line 205 "zsytrs_3.f"
    b_dim1 = *ldb;
#line 205 "zsytrs_3.f"
    b_offset = 1 + b_dim1;
#line 205 "zsytrs_3.f"
    b -= b_offset;
#line 205 "zsytrs_3.f"

#line 205 "zsytrs_3.f"
    /* Function Body */
#line 205 "zsytrs_3.f"
    *info = 0;
#line 206 "zsytrs_3.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 207 "zsytrs_3.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 208 "zsytrs_3.f"
	*info = -1;
#line 209 "zsytrs_3.f"
    } else if (*n < 0) {
#line 210 "zsytrs_3.f"
	*info = -2;
#line 211 "zsytrs_3.f"
    } else if (*nrhs < 0) {
#line 212 "zsytrs_3.f"
	*info = -3;
#line 213 "zsytrs_3.f"
    } else if (*lda < max(1,*n)) {
#line 214 "zsytrs_3.f"
	*info = -5;
#line 215 "zsytrs_3.f"
    } else if (*ldb < max(1,*n)) {
#line 216 "zsytrs_3.f"
	*info = -9;
#line 217 "zsytrs_3.f"
    }
#line 218 "zsytrs_3.f"
    if (*info != 0) {
#line 219 "zsytrs_3.f"
	i__1 = -(*info);
#line 219 "zsytrs_3.f"
	xerbla_("ZSYTRS_3", &i__1, (ftnlen)8);
#line 220 "zsytrs_3.f"
	return 0;
#line 221 "zsytrs_3.f"
    }

/*     Quick return if possible */

#line 225 "zsytrs_3.f"
    if (*n == 0 || *nrhs == 0) {
#line 225 "zsytrs_3.f"
	return 0;
#line 225 "zsytrs_3.f"
    }

#line 228 "zsytrs_3.f"
    if (upper) {

/*        Begin Upper */

/*        Solve A*X = B, where A = U*D*U**T. */

/*        P**T * B */

/*        Interchange rows K and IPIV(K) of matrix B in the same order */
/*        that the formation order of IPIV(I) vector for Upper case. */

/*        (We can do the simple loop over IPIV with decrement -1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 243 "zsytrs_3.f"
	for (k = *n; k >= 1; --k) {
#line 244 "zsytrs_3.f"
	    kp = (i__1 = ipiv[k], abs(i__1));
#line 245 "zsytrs_3.f"
	    if (kp != k) {
#line 246 "zsytrs_3.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 247 "zsytrs_3.f"
	    }
#line 248 "zsytrs_3.f"
	}

/*        Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 252 "zsytrs_3.f"
	ztrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 256 "zsytrs_3.f"
	i__ = *n;
#line 257 "zsytrs_3.f"
	while(i__ >= 1) {
#line 258 "zsytrs_3.f"
	    if (ipiv[i__] > 0) {
#line 259 "zsytrs_3.f"
		z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
#line 259 "zsytrs_3.f"
		zscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
#line 260 "zsytrs_3.f"
	    } else if (i__ > 1) {
#line 261 "zsytrs_3.f"
		i__1 = i__;
#line 261 "zsytrs_3.f"
		akm1k.r = e[i__1].r, akm1k.i = e[i__1].i;
#line 262 "zsytrs_3.f"
		z_div(&z__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
#line 262 "zsytrs_3.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 263 "zsytrs_3.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
#line 263 "zsytrs_3.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 264 "zsytrs_3.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 264 "zsytrs_3.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 264 "zsytrs_3.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 265 "zsytrs_3.f"
		i__1 = *nrhs;
#line 265 "zsytrs_3.f"
		for (j = 1; j <= i__1; ++j) {
#line 266 "zsytrs_3.f"
		    z_div(&z__1, &b[i__ - 1 + j * b_dim1], &akm1k);
#line 266 "zsytrs_3.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 267 "zsytrs_3.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
#line 267 "zsytrs_3.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 268 "zsytrs_3.f"
		    i__2 = i__ - 1 + j * b_dim1;
#line 268 "zsytrs_3.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 268 "zsytrs_3.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 268 "zsytrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 268 "zsytrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 269 "zsytrs_3.f"
		    i__2 = i__ + j * b_dim1;
#line 269 "zsytrs_3.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 269 "zsytrs_3.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 269 "zsytrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 269 "zsytrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 270 "zsytrs_3.f"
		}
#line 271 "zsytrs_3.f"
		--i__;
#line 272 "zsytrs_3.f"
	    }
#line 273 "zsytrs_3.f"
	    --i__;
#line 274 "zsytrs_3.f"
	}

/*        Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ] */

#line 278 "zsytrs_3.f"
	ztrsm_("L", "U", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ] */

/*        Interchange rows K and IPIV(K) of matrix B in reverse order */
/*        from the formation order of IPIV(I) vector for Upper case. */

/*        (We can do the simple loop over IPIV with increment 1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 289 "zsytrs_3.f"
	i__1 = *n;
#line 289 "zsytrs_3.f"
	for (k = 1; k <= i__1; ++k) {
#line 290 "zsytrs_3.f"
	    kp = (i__2 = ipiv[k], abs(i__2));
#line 291 "zsytrs_3.f"
	    if (kp != k) {
#line 292 "zsytrs_3.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 293 "zsytrs_3.f"
	    }
#line 294 "zsytrs_3.f"
	}

#line 296 "zsytrs_3.f"
    } else {

/*        Begin Lower */

/*        Solve A*X = B, where A = L*D*L**T. */

/*        P**T * B */
/*        Interchange rows K and IPIV(K) of matrix B in the same order */
/*        that the formation order of IPIV(I) vector for Lower case. */

/*        (We can do the simple loop over IPIV with increment 1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 310 "zsytrs_3.f"
	i__1 = *n;
#line 310 "zsytrs_3.f"
	for (k = 1; k <= i__1; ++k) {
#line 311 "zsytrs_3.f"
	    kp = (i__2 = ipiv[k], abs(i__2));
#line 312 "zsytrs_3.f"
	    if (kp != k) {
#line 313 "zsytrs_3.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 314 "zsytrs_3.f"
	    }
#line 315 "zsytrs_3.f"
	}

/*        Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 319 "zsytrs_3.f"
	ztrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 323 "zsytrs_3.f"
	i__ = 1;
#line 324 "zsytrs_3.f"
	while(i__ <= *n) {
#line 325 "zsytrs_3.f"
	    if (ipiv[i__] > 0) {
#line 326 "zsytrs_3.f"
		z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
#line 326 "zsytrs_3.f"
		zscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
#line 327 "zsytrs_3.f"
	    } else if (i__ < *n) {
#line 328 "zsytrs_3.f"
		i__1 = i__;
#line 328 "zsytrs_3.f"
		akm1k.r = e[i__1].r, akm1k.i = e[i__1].i;
#line 329 "zsytrs_3.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
#line 329 "zsytrs_3.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 330 "zsytrs_3.f"
		z_div(&z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
#line 330 "zsytrs_3.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 331 "zsytrs_3.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 331 "zsytrs_3.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 331 "zsytrs_3.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 332 "zsytrs_3.f"
		i__1 = *nrhs;
#line 332 "zsytrs_3.f"
		for (j = 1; j <= i__1; ++j) {
#line 333 "zsytrs_3.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
#line 333 "zsytrs_3.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 334 "zsytrs_3.f"
		    z_div(&z__1, &b[i__ + 1 + j * b_dim1], &akm1k);
#line 334 "zsytrs_3.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 335 "zsytrs_3.f"
		    i__2 = i__ + j * b_dim1;
#line 335 "zsytrs_3.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 335 "zsytrs_3.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 335 "zsytrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 335 "zsytrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 336 "zsytrs_3.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 336 "zsytrs_3.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 336 "zsytrs_3.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 336 "zsytrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 336 "zsytrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 337 "zsytrs_3.f"
		}
#line 338 "zsytrs_3.f"
		++i__;
#line 339 "zsytrs_3.f"
	    }
#line 340 "zsytrs_3.f"
	    ++i__;
#line 341 "zsytrs_3.f"
	}

/*        Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ] */

#line 345 "zsytrs_3.f"
	ztrsm_("L", "L", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ] */

/*        Interchange rows K and IPIV(K) of matrix B in reverse order */
/*        from the formation order of IPIV(I) vector for Lower case. */

/*        (We can do the simple loop over IPIV with decrement -1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 356 "zsytrs_3.f"
	for (k = *n; k >= 1; --k) {
#line 357 "zsytrs_3.f"
	    kp = (i__1 = ipiv[k], abs(i__1));
#line 358 "zsytrs_3.f"
	    if (kp != k) {
#line 359 "zsytrs_3.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 360 "zsytrs_3.f"
	    }
#line 361 "zsytrs_3.f"
	}

/*        END Lower */

#line 365 "zsytrs_3.f"
    }

#line 367 "zsytrs_3.f"
    return 0;

/*     End of ZSYTRS_3 */

} /* zsytrs_3__ */

