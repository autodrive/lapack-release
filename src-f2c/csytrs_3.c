#line 1 "csytrs_3.f"
/* csytrs_3.f -- translated by f2c (version 20100827).
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

#line 1 "csytrs_3.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b CSYTRS_3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRS_3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs_
3.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs_
3.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs_
3.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, */
/*                            INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > CSYTRS_3 solves a system of linear equations A * X = B with a complex */
/* > symmetric matrix A using the factorization computed */
/* > by CSYTRF_RK or CSYTRF_BK: */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          Diagonal of the block diagonal matrix D and factors U or L */
/* >          as computed by CSYTRF_RK and CSYTRF_BK: */
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
/* >          E is COMPLEX array, dimension (N) */
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
/* >          as determined by CSYTRF_RK or CSYTRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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

/* > \ingroup complexSYcomputational */

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
/* Subroutine */ int csytrs_3__(char *uplo, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex denom;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ctrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

#line 205 "csytrs_3.f"
    /* Parameter adjustments */
#line 205 "csytrs_3.f"
    a_dim1 = *lda;
#line 205 "csytrs_3.f"
    a_offset = 1 + a_dim1;
#line 205 "csytrs_3.f"
    a -= a_offset;
#line 205 "csytrs_3.f"
    --e;
#line 205 "csytrs_3.f"
    --ipiv;
#line 205 "csytrs_3.f"
    b_dim1 = *ldb;
#line 205 "csytrs_3.f"
    b_offset = 1 + b_dim1;
#line 205 "csytrs_3.f"
    b -= b_offset;
#line 205 "csytrs_3.f"

#line 205 "csytrs_3.f"
    /* Function Body */
#line 205 "csytrs_3.f"
    *info = 0;
#line 206 "csytrs_3.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 207 "csytrs_3.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 208 "csytrs_3.f"
	*info = -1;
#line 209 "csytrs_3.f"
    } else if (*n < 0) {
#line 210 "csytrs_3.f"
	*info = -2;
#line 211 "csytrs_3.f"
    } else if (*nrhs < 0) {
#line 212 "csytrs_3.f"
	*info = -3;
#line 213 "csytrs_3.f"
    } else if (*lda < max(1,*n)) {
#line 214 "csytrs_3.f"
	*info = -5;
#line 215 "csytrs_3.f"
    } else if (*ldb < max(1,*n)) {
#line 216 "csytrs_3.f"
	*info = -9;
#line 217 "csytrs_3.f"
    }
#line 218 "csytrs_3.f"
    if (*info != 0) {
#line 219 "csytrs_3.f"
	i__1 = -(*info);
#line 219 "csytrs_3.f"
	xerbla_("CSYTRS_3", &i__1, (ftnlen)8);
#line 220 "csytrs_3.f"
	return 0;
#line 221 "csytrs_3.f"
    }

/*     Quick return if possible */

#line 225 "csytrs_3.f"
    if (*n == 0 || *nrhs == 0) {
#line 225 "csytrs_3.f"
	return 0;
#line 225 "csytrs_3.f"
    }

#line 228 "csytrs_3.f"
    if (upper) {

/*        Begin Upper */

/*        Solve A*X = B, where A = U*D*U**T. */

/*        P**T * B */

/*        Interchange rows K and IPIV(K) of matrix B in the same order */
/*        that the formation order of IPIV(I) vector for Upper case. */

/*        (We can do the simple loop over IPIV with decrement -1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 243 "csytrs_3.f"
	for (k = *n; k >= 1; --k) {
#line 244 "csytrs_3.f"
	    kp = (i__1 = ipiv[k], abs(i__1));
#line 245 "csytrs_3.f"
	    if (kp != k) {
#line 246 "csytrs_3.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 247 "csytrs_3.f"
	    }
#line 248 "csytrs_3.f"
	}

/*        Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 252 "csytrs_3.f"
	ctrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 256 "csytrs_3.f"
	i__ = *n;
#line 257 "csytrs_3.f"
	while(i__ >= 1) {
#line 258 "csytrs_3.f"
	    if (ipiv[i__] > 0) {
#line 259 "csytrs_3.f"
		z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
#line 259 "csytrs_3.f"
		cscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
#line 260 "csytrs_3.f"
	    } else if (i__ > 1) {
#line 261 "csytrs_3.f"
		i__1 = i__;
#line 261 "csytrs_3.f"
		akm1k.r = e[i__1].r, akm1k.i = e[i__1].i;
#line 262 "csytrs_3.f"
		z_div(&z__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
#line 262 "csytrs_3.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 263 "csytrs_3.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
#line 263 "csytrs_3.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 264 "csytrs_3.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 264 "csytrs_3.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 264 "csytrs_3.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 265 "csytrs_3.f"
		i__1 = *nrhs;
#line 265 "csytrs_3.f"
		for (j = 1; j <= i__1; ++j) {
#line 266 "csytrs_3.f"
		    z_div(&z__1, &b[i__ - 1 + j * b_dim1], &akm1k);
#line 266 "csytrs_3.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 267 "csytrs_3.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
#line 267 "csytrs_3.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 268 "csytrs_3.f"
		    i__2 = i__ - 1 + j * b_dim1;
#line 268 "csytrs_3.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 268 "csytrs_3.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 268 "csytrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 268 "csytrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 269 "csytrs_3.f"
		    i__2 = i__ + j * b_dim1;
#line 269 "csytrs_3.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 269 "csytrs_3.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 269 "csytrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 269 "csytrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 270 "csytrs_3.f"
		}
#line 271 "csytrs_3.f"
		--i__;
#line 272 "csytrs_3.f"
	    }
#line 273 "csytrs_3.f"
	    --i__;
#line 274 "csytrs_3.f"
	}

/*        Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ] */

#line 278 "csytrs_3.f"
	ctrsm_("L", "U", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ] */

/*        Interchange rows K and IPIV(K) of matrix B in reverse order */
/*        from the formation order of IPIV(I) vector for Upper case. */

/*        (We can do the simple loop over IPIV with increment 1, */
/*        since the ABS value of IPIV( I ) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 289 "csytrs_3.f"
	i__1 = *n;
#line 289 "csytrs_3.f"
	for (k = 1; k <= i__1; ++k) {
#line 290 "csytrs_3.f"
	    kp = (i__2 = ipiv[k], abs(i__2));
#line 291 "csytrs_3.f"
	    if (kp != k) {
#line 292 "csytrs_3.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 293 "csytrs_3.f"
	    }
#line 294 "csytrs_3.f"
	}

#line 296 "csytrs_3.f"
    } else {

/*        Begin Lower */

/*        Solve A*X = B, where A = L*D*L**T. */

/*        P**T * B */
/*        Interchange rows K and IPIV(K) of matrix B in the same order */
/*        that the formation order of IPIV(I) vector for Lower case. */

/*        (We can do the simple loop over IPIV with increment 1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 310 "csytrs_3.f"
	i__1 = *n;
#line 310 "csytrs_3.f"
	for (k = 1; k <= i__1; ++k) {
#line 311 "csytrs_3.f"
	    kp = (i__2 = ipiv[k], abs(i__2));
#line 312 "csytrs_3.f"
	    if (kp != k) {
#line 313 "csytrs_3.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 314 "csytrs_3.f"
	    }
#line 315 "csytrs_3.f"
	}

/*        Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 319 "csytrs_3.f"
	ctrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 323 "csytrs_3.f"
	i__ = 1;
#line 324 "csytrs_3.f"
	while(i__ <= *n) {
#line 325 "csytrs_3.f"
	    if (ipiv[i__] > 0) {
#line 326 "csytrs_3.f"
		z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
#line 326 "csytrs_3.f"
		cscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
#line 327 "csytrs_3.f"
	    } else if (i__ < *n) {
#line 328 "csytrs_3.f"
		i__1 = i__;
#line 328 "csytrs_3.f"
		akm1k.r = e[i__1].r, akm1k.i = e[i__1].i;
#line 329 "csytrs_3.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
#line 329 "csytrs_3.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 330 "csytrs_3.f"
		z_div(&z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
#line 330 "csytrs_3.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 331 "csytrs_3.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 331 "csytrs_3.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 331 "csytrs_3.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 332 "csytrs_3.f"
		i__1 = *nrhs;
#line 332 "csytrs_3.f"
		for (j = 1; j <= i__1; ++j) {
#line 333 "csytrs_3.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
#line 333 "csytrs_3.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 334 "csytrs_3.f"
		    z_div(&z__1, &b[i__ + 1 + j * b_dim1], &akm1k);
#line 334 "csytrs_3.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 335 "csytrs_3.f"
		    i__2 = i__ + j * b_dim1;
#line 335 "csytrs_3.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 335 "csytrs_3.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 335 "csytrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 335 "csytrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 336 "csytrs_3.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 336 "csytrs_3.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 336 "csytrs_3.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 336 "csytrs_3.f"
		    z_div(&z__1, &z__2, &denom);
#line 336 "csytrs_3.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 337 "csytrs_3.f"
		}
#line 338 "csytrs_3.f"
		++i__;
#line 339 "csytrs_3.f"
	    }
#line 340 "csytrs_3.f"
	    ++i__;
#line 341 "csytrs_3.f"
	}

/*        Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ] */

#line 345 "csytrs_3.f"
	ctrsm_("L", "L", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ] */

/*        Interchange rows K and IPIV(K) of matrix B in reverse order */
/*        from the formation order of IPIV(I) vector for Lower case. */

/*        (We can do the simple loop over IPIV with decrement -1, */
/*        since the ABS value of IPIV(I) represents the row index */
/*        of the interchange with row i in both 1x1 and 2x2 pivot cases) */

#line 356 "csytrs_3.f"
	for (k = *n; k >= 1; --k) {
#line 357 "csytrs_3.f"
	    kp = (i__1 = ipiv[k], abs(i__1));
#line 358 "csytrs_3.f"
	    if (kp != k) {
#line 359 "csytrs_3.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 360 "csytrs_3.f"
	    }
#line 361 "csytrs_3.f"
	}

/*        END Lower */

#line 365 "csytrs_3.f"
    }

#line 367 "csytrs_3.f"
    return 0;

/*     End of CSYTRS_3 */

} /* csytrs_3__ */

