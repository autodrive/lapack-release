#line 1 "ssytri_3x.f"
/* ssytri_3x.f -- translated by f2c (version 20100827).
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

#line 1 "ssytri_3x.f"
/* Table of constant values */

static doublereal c_b10 = 1.;
static doublereal c_b14 = 0.;

/* > \brief \b SSYTRI_3X */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTRI_3X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri_
3x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri_
3x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri_
3x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ),  E( * ), WORK( N+NB+1, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > SSYTRI_3X computes the inverse of a real symmetric indefinite */
/* > matrix A using the factorization computed by SSYTRF_RK or SSYTRF_BK: */
/* > */
/* >     A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is symmetric and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix. */
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
/* >          On entry, diagonal of the block diagonal matrix D and */
/* >          factors U or L as computed by SYTRF_RK and SSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                should be provided on entry in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, if INFO = 0, the symmetric inverse of the original */
/* >          matrix. */
/* >             If UPLO = 'U': the upper triangular part of the inverse */
/* >             is formed and the part of A below the diagonal is not */
/* >             referenced; */
/* >             If UPLO = 'L': the lower triangular part of the inverse */
/* >             is formed and the part of A above the diagonal is not */
/* >             referenced. */
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
/* >          E is REAL array, dimension (N) */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) not referenced. */
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
/* >          as determined by SSYTRF_RK or SSYTRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N+NB+1,NB+3). */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          Block size. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its */
/* >               inverse could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup singleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > \verbatim */
/* > */
/* >  June 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int ssytri_3x__(char *uplo, integer *n, doublereal *a, 
	integer *lda, doublereal *e, integer *ipiv, doublereal *work, integer 
	*nb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, k;
    extern /* Subroutine */ int ssyswapr_(char *, integer *, doublereal *, 
	    integer *, integer *, integer *, ftnlen);
    static doublereal t, ak;
    static integer u11, ip, nnb, cut;
    static doublereal akp1;
    static integer invd;
    static doublereal akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal u01_i_j__, u11_i_j__;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer icount;
    extern /* Subroutine */ int strtri_(char *, char *, integer *, doublereal 
	    *, integer *, integer *, ftnlen, ftnlen);
    static doublereal u01_ip1_j__, u11_ip1_j__;


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

/*     Test the input parameters. */

#line 202 "ssytri_3x.f"
    /* Parameter adjustments */
#line 202 "ssytri_3x.f"
    a_dim1 = *lda;
#line 202 "ssytri_3x.f"
    a_offset = 1 + a_dim1;
#line 202 "ssytri_3x.f"
    a -= a_offset;
#line 202 "ssytri_3x.f"
    --e;
#line 202 "ssytri_3x.f"
    --ipiv;
#line 202 "ssytri_3x.f"
    work_dim1 = *n + *nb + 1;
#line 202 "ssytri_3x.f"
    work_offset = 1 + work_dim1;
#line 202 "ssytri_3x.f"
    work -= work_offset;
#line 202 "ssytri_3x.f"

#line 202 "ssytri_3x.f"
    /* Function Body */
#line 202 "ssytri_3x.f"
    *info = 0;
#line 203 "ssytri_3x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 204 "ssytri_3x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 205 "ssytri_3x.f"
	*info = -1;
#line 206 "ssytri_3x.f"
    } else if (*n < 0) {
#line 207 "ssytri_3x.f"
	*info = -2;
#line 208 "ssytri_3x.f"
    } else if (*lda < max(1,*n)) {
#line 209 "ssytri_3x.f"
	*info = -4;
#line 210 "ssytri_3x.f"
    }

/*     Quick return if possible */

#line 214 "ssytri_3x.f"
    if (*info != 0) {
#line 215 "ssytri_3x.f"
	i__1 = -(*info);
#line 215 "ssytri_3x.f"
	xerbla_("SSYTRI_3X", &i__1, (ftnlen)9);
#line 216 "ssytri_3x.f"
	return 0;
#line 217 "ssytri_3x.f"
    }
#line 218 "ssytri_3x.f"
    if (*n == 0) {
#line 218 "ssytri_3x.f"
	return 0;
#line 218 "ssytri_3x.f"
    }

/*     Workspace got Non-diag elements of D */

#line 223 "ssytri_3x.f"
    i__1 = *n;
#line 223 "ssytri_3x.f"
    for (k = 1; k <= i__1; ++k) {
#line 224 "ssytri_3x.f"
	work[k + work_dim1] = e[k];
#line 225 "ssytri_3x.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 229 "ssytri_3x.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 233 "ssytri_3x.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 234 "ssytri_3x.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 234 "ssytri_3x.f"
		return 0;
#line 234 "ssytri_3x.f"
	    }
#line 236 "ssytri_3x.f"
	}
#line 237 "ssytri_3x.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 241 "ssytri_3x.f"
	i__1 = *n;
#line 241 "ssytri_3x.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 242 "ssytri_3x.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 242 "ssytri_3x.f"
		return 0;
#line 242 "ssytri_3x.f"
	    }
#line 244 "ssytri_3x.f"
	}
#line 245 "ssytri_3x.f"
    }

#line 247 "ssytri_3x.f"
    *info = 0;

/*     Splitting Workspace */
/*     U01 is a block ( N, NB+1 ) */
/*     The first element of U01 is in WORK( 1, 1 ) */
/*     U11 is a block ( NB+1, NB+1 ) */
/*     The first element of U11 is in WORK( N+1, 1 ) */

#line 255 "ssytri_3x.f"
    u11 = *n;

/*     INVD is a block ( N, 2 ) */
/*     The first element of INVD is in WORK( 1, INVD ) */

#line 260 "ssytri_3x.f"
    invd = *nb + 2;
#line 262 "ssytri_3x.f"
    if (upper) {

/*        Begin Upper */

/*        invA = P * inv(U**T) * inv(D) * inv(U) * P**T. */

#line 268 "ssytri_3x.f"
	strtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*        inv(D) and inv(D) * inv(U) */

#line 272 "ssytri_3x.f"
	k = 1;
#line 273 "ssytri_3x.f"
	while(k <= *n) {
#line 274 "ssytri_3x.f"
	    if (ipiv[k] > 0) {
/*              1 x 1 diagonal NNB */
#line 276 "ssytri_3x.f"
		work[k + invd * work_dim1] = 1. / a[k + k * a_dim1];
#line 277 "ssytri_3x.f"
		work[k + (invd + 1) * work_dim1] = 0.;
#line 278 "ssytri_3x.f"
	    } else {
/*              2 x 2 diagonal NNB */
#line 280 "ssytri_3x.f"
		t = work[k + 1 + work_dim1];
#line 281 "ssytri_3x.f"
		ak = a[k + k * a_dim1] / t;
#line 282 "ssytri_3x.f"
		akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
#line 283 "ssytri_3x.f"
		akkp1 = work[k + 1 + work_dim1] / t;
#line 284 "ssytri_3x.f"
		d__ = t * (ak * akp1 - 1.);
#line 285 "ssytri_3x.f"
		work[k + invd * work_dim1] = akp1 / d__;
#line 286 "ssytri_3x.f"
		work[k + 1 + (invd + 1) * work_dim1] = ak / d__;
#line 287 "ssytri_3x.f"
		work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
#line 288 "ssytri_3x.f"
		work[k + 1 + invd * work_dim1] = work[k + (invd + 1) * 
			work_dim1];
#line 289 "ssytri_3x.f"
		++k;
#line 290 "ssytri_3x.f"
	    }
#line 291 "ssytri_3x.f"
	    ++k;
#line 292 "ssytri_3x.f"
	}

/*        inv(U**T) = (inv(U))**T */

/*        inv(U**T) * inv(D) * inv(U) */

#line 298 "ssytri_3x.f"
	cut = *n;
#line 299 "ssytri_3x.f"
	while(cut > 0) {
#line 300 "ssytri_3x.f"
	    nnb = *nb;
#line 301 "ssytri_3x.f"
	    if (cut <= nnb) {
#line 302 "ssytri_3x.f"
		nnb = cut;
#line 303 "ssytri_3x.f"
	    } else {
#line 304 "ssytri_3x.f"
		icount = 0;
/*              count negative elements, */
#line 306 "ssytri_3x.f"
		i__1 = cut;
#line 306 "ssytri_3x.f"
		for (i__ = cut + 1 - nnb; i__ <= i__1; ++i__) {
#line 307 "ssytri_3x.f"
		    if (ipiv[i__] < 0) {
#line 307 "ssytri_3x.f"
			++icount;
#line 307 "ssytri_3x.f"
		    }
#line 308 "ssytri_3x.f"
		}
/*              need a even number for a clear cut */
#line 310 "ssytri_3x.f"
		if (icount % 2 == 1) {
#line 310 "ssytri_3x.f"
		    ++nnb;
#line 310 "ssytri_3x.f"
		}
#line 311 "ssytri_3x.f"
	    }
#line 313 "ssytri_3x.f"
	    cut -= nnb;

/*           U01 Block */

#line 317 "ssytri_3x.f"
	    i__1 = cut;
#line 317 "ssytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 318 "ssytri_3x.f"
		i__2 = nnb;
#line 318 "ssytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 319 "ssytri_3x.f"
		    work[i__ + j * work_dim1] = a[i__ + (cut + j) * a_dim1];
#line 320 "ssytri_3x.f"
		}
#line 321 "ssytri_3x.f"
	    }

/*           U11 Block */

#line 325 "ssytri_3x.f"
	    i__1 = nnb;
#line 325 "ssytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 326 "ssytri_3x.f"
		work[u11 + i__ + i__ * work_dim1] = 1.;
#line 327 "ssytri_3x.f"
		i__2 = i__ - 1;
#line 327 "ssytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 328 "ssytri_3x.f"
		    work[u11 + i__ + j * work_dim1] = 0.;
#line 329 "ssytri_3x.f"
		}
#line 330 "ssytri_3x.f"
		i__2 = nnb;
#line 330 "ssytri_3x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 331 "ssytri_3x.f"
		    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) 
			    * a_dim1];
#line 332 "ssytri_3x.f"
		}
#line 333 "ssytri_3x.f"
	    }

/*           invD * U01 */

#line 337 "ssytri_3x.f"
	    i__ = 1;
#line 338 "ssytri_3x.f"
	    while(i__ <= cut) {
#line 339 "ssytri_3x.f"
		if (ipiv[i__] > 0) {
#line 340 "ssytri_3x.f"
		    i__1 = nnb;
#line 340 "ssytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 341 "ssytri_3x.f"
			work[i__ + j * work_dim1] = work[i__ + invd * 
				work_dim1] * work[i__ + j * work_dim1];
#line 342 "ssytri_3x.f"
		    }
#line 343 "ssytri_3x.f"
		} else {
#line 344 "ssytri_3x.f"
		    i__1 = nnb;
#line 344 "ssytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 345 "ssytri_3x.f"
			u01_i_j__ = work[i__ + j * work_dim1];
#line 346 "ssytri_3x.f"
			u01_ip1_j__ = work[i__ + 1 + j * work_dim1];
#line 347 "ssytri_3x.f"
			work[i__ + j * work_dim1] = work[i__ + invd * 
				work_dim1] * u01_i_j__ + work[i__ + (invd + 1)
				 * work_dim1] * u01_ip1_j__;
#line 349 "ssytri_3x.f"
			work[i__ + 1 + j * work_dim1] = work[i__ + 1 + invd * 
				work_dim1] * u01_i_j__ + work[i__ + 1 + (invd 
				+ 1) * work_dim1] * u01_ip1_j__;
#line 351 "ssytri_3x.f"
		    }
#line 352 "ssytri_3x.f"
		    ++i__;
#line 353 "ssytri_3x.f"
		}
#line 354 "ssytri_3x.f"
		++i__;
#line 355 "ssytri_3x.f"
	    }

/*           invD1 * U11 */

#line 359 "ssytri_3x.f"
	    i__ = 1;
#line 360 "ssytri_3x.f"
	    while(i__ <= nnb) {
#line 361 "ssytri_3x.f"
		if (ipiv[cut + i__] > 0) {
#line 362 "ssytri_3x.f"
		    i__1 = nnb;
#line 362 "ssytri_3x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 363 "ssytri_3x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1];
#line 364 "ssytri_3x.f"
		    }
#line 365 "ssytri_3x.f"
		} else {
#line 366 "ssytri_3x.f"
		    i__1 = nnb;
#line 366 "ssytri_3x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 367 "ssytri_3x.f"
			u11_i_j__ = work[u11 + i__ + j * work_dim1];
#line 368 "ssytri_3x.f"
			u11_ip1_j__ = work[u11 + i__ + 1 + j * work_dim1];
#line 369 "ssytri_3x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1] + work[cut + i__ + (invd + 1) * 
				work_dim1] * work[u11 + i__ + 1 + j * 
				work_dim1];
#line 371 "ssytri_3x.f"
			work[u11 + i__ + 1 + j * work_dim1] = work[cut + i__ 
				+ 1 + invd * work_dim1] * u11_i_j__ + work[
				cut + i__ + 1 + (invd + 1) * work_dim1] * 
				u11_ip1_j__;
#line 373 "ssytri_3x.f"
		    }
#line 374 "ssytri_3x.f"
		    ++i__;
#line 375 "ssytri_3x.f"
		}
#line 376 "ssytri_3x.f"
		++i__;
#line 377 "ssytri_3x.f"
	    }

/*           U11**T * invD1 * U11 -> U11 */

#line 381 "ssytri_3x.f"
	    i__1 = *n + *nb + 1;
#line 381 "ssytri_3x.f"
	    strmm_("L", "U", "T", "U", &nnb, &nnb, &c_b10, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 385 "ssytri_3x.f"
	    i__1 = nnb;
#line 385 "ssytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "ssytri_3x.f"
		i__2 = nnb;
#line 386 "ssytri_3x.f"
		for (j = i__; j <= i__2; ++j) {
#line 387 "ssytri_3x.f"
		    a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * 
			    work_dim1];
#line 388 "ssytri_3x.f"
		}
#line 389 "ssytri_3x.f"
	    }

/*           U01**T * invD * U01 -> A( CUT+I, CUT+J ) */

#line 393 "ssytri_3x.f"
	    i__1 = *n + *nb + 1;
#line 393 "ssytri_3x.f"
	    i__2 = *n + *nb + 1;
#line 393 "ssytri_3x.f"
	    sgemm_("T", "N", &nnb, &nnb, &cut, &c_b10, &a[(cut + 1) * a_dim1 
		    + 1], lda, &work[work_offset], &i__1, &c_b14, &work[u11 + 
		    1 + work_dim1], &i__2, (ftnlen)1, (ftnlen)1);

/*           U11 =  U11**T * invD1 * U11 + U01**T * invD * U01 */

#line 399 "ssytri_3x.f"
	    i__1 = nnb;
#line 399 "ssytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 400 "ssytri_3x.f"
		i__2 = nnb;
#line 400 "ssytri_3x.f"
		for (j = i__; j <= i__2; ++j) {
#line 401 "ssytri_3x.f"
		    a[cut + i__ + (cut + j) * a_dim1] += work[u11 + i__ + j * 
			    work_dim1];
#line 402 "ssytri_3x.f"
		}
#line 403 "ssytri_3x.f"
	    }

/*           U01 =  U00**T * invD0 * U01 */

#line 407 "ssytri_3x.f"
	    i__1 = *n + *nb + 1;
#line 407 "ssytri_3x.f"
	    strmm_("L", uplo, "T", "U", &cut, &nnb, &c_b10, &a[a_offset], lda,
		     &work[work_offset], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);

/*           Update U01 */

#line 413 "ssytri_3x.f"
	    i__1 = cut;
#line 413 "ssytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 414 "ssytri_3x.f"
		i__2 = nnb;
#line 414 "ssytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 415 "ssytri_3x.f"
		    a[i__ + (cut + j) * a_dim1] = work[i__ + j * work_dim1];
#line 416 "ssytri_3x.f"
		}
#line 417 "ssytri_3x.f"
	    }

/*           Next Block */

#line 421 "ssytri_3x.f"
	}

/*        Apply PERMUTATIONS P and P**T: */
/*        P * inv(U**T) * inv(D) * inv(U) * P**T. */
/*        Interchange rows and columns I and IPIV(I) in reverse order */
/*        from the formation order of IPIV vector for Upper case. */

/*        ( We can use a loop over IPIV with increment 1, */
/*        since the ABS value of IPIV(I) represents the row (column) */
/*        index of the interchange with row (column) i in both 1x1 */
/*        and 2x2 pivot cases, i.e. we don't need separate code branches */
/*        for 1x1 and 2x2 pivot cases ) */

#line 434 "ssytri_3x.f"
	i__1 = *n;
#line 434 "ssytri_3x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 435 "ssytri_3x.f"
	    ip = (i__2 = ipiv[i__], abs(i__2));
#line 436 "ssytri_3x.f"
	    if (ip != i__) {
#line 437 "ssytri_3x.f"
		if (i__ < ip) {
#line 437 "ssytri_3x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 437 "ssytri_3x.f"
		}
#line 438 "ssytri_3x.f"
		if (i__ > ip) {
#line 438 "ssytri_3x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 438 "ssytri_3x.f"
		}
#line 439 "ssytri_3x.f"
	    }
#line 440 "ssytri_3x.f"
	}

#line 442 "ssytri_3x.f"
    } else {

/*        Begin Lower */

/*        inv A = P * inv(L**T) * inv(D) * inv(L) * P**T. */

#line 448 "ssytri_3x.f"
	strtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*        inv(D) and inv(D) * inv(L) */

#line 452 "ssytri_3x.f"
	k = *n;
#line 453 "ssytri_3x.f"
	while(k >= 1) {
#line 454 "ssytri_3x.f"
	    if (ipiv[k] > 0) {
/*              1 x 1 diagonal NNB */
#line 456 "ssytri_3x.f"
		work[k + invd * work_dim1] = 1. / a[k + k * a_dim1];
#line 457 "ssytri_3x.f"
		work[k + (invd + 1) * work_dim1] = 0.;
#line 458 "ssytri_3x.f"
	    } else {
/*              2 x 2 diagonal NNB */
#line 460 "ssytri_3x.f"
		t = work[k - 1 + work_dim1];
#line 461 "ssytri_3x.f"
		ak = a[k - 1 + (k - 1) * a_dim1] / t;
#line 462 "ssytri_3x.f"
		akp1 = a[k + k * a_dim1] / t;
#line 463 "ssytri_3x.f"
		akkp1 = work[k - 1 + work_dim1] / t;
#line 464 "ssytri_3x.f"
		d__ = t * (ak * akp1 - 1.);
#line 465 "ssytri_3x.f"
		work[k - 1 + invd * work_dim1] = akp1 / d__;
#line 466 "ssytri_3x.f"
		work[k + invd * work_dim1] = ak / d__;
#line 467 "ssytri_3x.f"
		work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
#line 468 "ssytri_3x.f"
		work[k - 1 + (invd + 1) * work_dim1] = work[k + (invd + 1) * 
			work_dim1];
#line 469 "ssytri_3x.f"
		--k;
#line 470 "ssytri_3x.f"
	    }
#line 471 "ssytri_3x.f"
	    --k;
#line 472 "ssytri_3x.f"
	}

/*        inv(L**T) = (inv(L))**T */

/*        inv(L**T) * inv(D) * inv(L) */

#line 478 "ssytri_3x.f"
	cut = 0;
#line 479 "ssytri_3x.f"
	while(cut < *n) {
#line 480 "ssytri_3x.f"
	    nnb = *nb;
#line 481 "ssytri_3x.f"
	    if (cut + nnb > *n) {
#line 482 "ssytri_3x.f"
		nnb = *n - cut;
#line 483 "ssytri_3x.f"
	    } else {
#line 484 "ssytri_3x.f"
		icount = 0;
/*              count negative elements, */
#line 486 "ssytri_3x.f"
		i__1 = cut + nnb;
#line 486 "ssytri_3x.f"
		for (i__ = cut + 1; i__ <= i__1; ++i__) {
#line 487 "ssytri_3x.f"
		    if (ipiv[i__] < 0) {
#line 487 "ssytri_3x.f"
			++icount;
#line 487 "ssytri_3x.f"
		    }
#line 488 "ssytri_3x.f"
		}
/*              need a even number for a clear cut */
#line 490 "ssytri_3x.f"
		if (icount % 2 == 1) {
#line 490 "ssytri_3x.f"
		    ++nnb;
#line 490 "ssytri_3x.f"
		}
#line 491 "ssytri_3x.f"
	    }

/*           L21 Block */

#line 495 "ssytri_3x.f"
	    i__1 = *n - cut - nnb;
#line 495 "ssytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 496 "ssytri_3x.f"
		i__2 = nnb;
#line 496 "ssytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 497 "ssytri_3x.f"
		    work[i__ + j * work_dim1] = a[cut + nnb + i__ + (cut + j) 
			    * a_dim1];
#line 498 "ssytri_3x.f"
		}
#line 499 "ssytri_3x.f"
	    }

/*           L11 Block */

#line 503 "ssytri_3x.f"
	    i__1 = nnb;
#line 503 "ssytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 504 "ssytri_3x.f"
		work[u11 + i__ + i__ * work_dim1] = 1.;
#line 505 "ssytri_3x.f"
		i__2 = nnb;
#line 505 "ssytri_3x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 506 "ssytri_3x.f"
		    work[u11 + i__ + j * work_dim1] = 0.;
#line 507 "ssytri_3x.f"
		}
#line 508 "ssytri_3x.f"
		i__2 = i__ - 1;
#line 508 "ssytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 509 "ssytri_3x.f"
		    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) 
			    * a_dim1];
#line 510 "ssytri_3x.f"
		}
#line 511 "ssytri_3x.f"
	    }

/*           invD*L21 */

#line 515 "ssytri_3x.f"
	    i__ = *n - cut - nnb;
#line 516 "ssytri_3x.f"
	    while(i__ >= 1) {
#line 517 "ssytri_3x.f"
		if (ipiv[cut + nnb + i__] > 0) {
#line 518 "ssytri_3x.f"
		    i__1 = nnb;
#line 518 "ssytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 519 "ssytri_3x.f"
			work[i__ + j * work_dim1] = work[cut + nnb + i__ + 
				invd * work_dim1] * work[i__ + j * work_dim1];
#line 520 "ssytri_3x.f"
		    }
#line 521 "ssytri_3x.f"
		} else {
#line 522 "ssytri_3x.f"
		    i__1 = nnb;
#line 522 "ssytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 523 "ssytri_3x.f"
			u01_i_j__ = work[i__ + j * work_dim1];
#line 524 "ssytri_3x.f"
			u01_ip1_j__ = work[i__ - 1 + j * work_dim1];
#line 525 "ssytri_3x.f"
			work[i__ + j * work_dim1] = work[cut + nnb + i__ + 
				invd * work_dim1] * u01_i_j__ + work[cut + 
				nnb + i__ + (invd + 1) * work_dim1] * 
				u01_ip1_j__;
#line 527 "ssytri_3x.f"
			work[i__ - 1 + j * work_dim1] = work[cut + nnb + i__ 
				- 1 + (invd + 1) * work_dim1] * u01_i_j__ + 
				work[cut + nnb + i__ - 1 + invd * work_dim1] *
				 u01_ip1_j__;
#line 529 "ssytri_3x.f"
		    }
#line 530 "ssytri_3x.f"
		    --i__;
#line 531 "ssytri_3x.f"
		}
#line 532 "ssytri_3x.f"
		--i__;
#line 533 "ssytri_3x.f"
	    }

/*           invD1*L11 */

#line 537 "ssytri_3x.f"
	    i__ = nnb;
#line 538 "ssytri_3x.f"
	    while(i__ >= 1) {
#line 539 "ssytri_3x.f"
		if (ipiv[cut + i__] > 0) {
#line 540 "ssytri_3x.f"
		    i__1 = nnb;
#line 540 "ssytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 541 "ssytri_3x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1];
#line 542 "ssytri_3x.f"
		    }
#line 544 "ssytri_3x.f"
		} else {
#line 545 "ssytri_3x.f"
		    i__1 = nnb;
#line 545 "ssytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 546 "ssytri_3x.f"
			u11_i_j__ = work[u11 + i__ + j * work_dim1];
#line 547 "ssytri_3x.f"
			u11_ip1_j__ = work[u11 + i__ - 1 + j * work_dim1];
#line 548 "ssytri_3x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1] + work[cut + i__ + (invd + 1) * 
				work_dim1] * u11_ip1_j__;
#line 550 "ssytri_3x.f"
			work[u11 + i__ - 1 + j * work_dim1] = work[cut + i__ 
				- 1 + (invd + 1) * work_dim1] * u11_i_j__ + 
				work[cut + i__ - 1 + invd * work_dim1] * 
				u11_ip1_j__;
#line 552 "ssytri_3x.f"
		    }
#line 553 "ssytri_3x.f"
		    --i__;
#line 554 "ssytri_3x.f"
		}
#line 555 "ssytri_3x.f"
		--i__;
#line 556 "ssytri_3x.f"
	    }

/*           L11**T * invD1 * L11 -> L11 */

#line 560 "ssytri_3x.f"
	    i__1 = *n + *nb + 1;
#line 560 "ssytri_3x.f"
	    strmm_("L", uplo, "T", "U", &nnb, &nnb, &c_b10, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 565 "ssytri_3x.f"
	    i__1 = nnb;
#line 565 "ssytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 566 "ssytri_3x.f"
		i__2 = i__;
#line 566 "ssytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 567 "ssytri_3x.f"
		    a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * 
			    work_dim1];
#line 568 "ssytri_3x.f"
		}
#line 569 "ssytri_3x.f"
	    }

#line 571 "ssytri_3x.f"
	    if (cut + nnb < *n) {

/*              L21**T * invD2*L21 -> A( CUT+I, CUT+J ) */

#line 575 "ssytri_3x.f"
		i__1 = *n - nnb - cut;
#line 575 "ssytri_3x.f"
		i__2 = *n + *nb + 1;
#line 575 "ssytri_3x.f"
		i__3 = *n + *nb + 1;
#line 575 "ssytri_3x.f"
		sgemm_("T", "N", &nnb, &nnb, &i__1, &c_b10, &a[cut + nnb + 1 
			+ (cut + 1) * a_dim1], lda, &work[work_offset], &i__2,
			 &c_b14, &work[u11 + 1 + work_dim1], &i__3, (ftnlen)1,
			 (ftnlen)1);

/*              L11 =  L11**T * invD1 * L11 + U01**T * invD * U01 */

#line 582 "ssytri_3x.f"
		i__1 = nnb;
#line 582 "ssytri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 583 "ssytri_3x.f"
		    i__2 = i__;
#line 583 "ssytri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 584 "ssytri_3x.f"
			a[cut + i__ + (cut + j) * a_dim1] += work[u11 + i__ + 
				j * work_dim1];
#line 585 "ssytri_3x.f"
		    }
#line 586 "ssytri_3x.f"
		}

/*              L01 =  L22**T * invD2 * L21 */

#line 590 "ssytri_3x.f"
		i__1 = *n - nnb - cut;
#line 590 "ssytri_3x.f"
		i__2 = *n + *nb + 1;
#line 590 "ssytri_3x.f"
		strmm_("L", uplo, "T", "U", &i__1, &nnb, &c_b10, &a[cut + nnb 
			+ 1 + (cut + nnb + 1) * a_dim1], lda, &work[
			work_offset], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);

/*              Update L21 */

#line 596 "ssytri_3x.f"
		i__1 = *n - cut - nnb;
#line 596 "ssytri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 597 "ssytri_3x.f"
		    i__2 = nnb;
#line 597 "ssytri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 598 "ssytri_3x.f"
			a[cut + nnb + i__ + (cut + j) * a_dim1] = work[i__ + 
				j * work_dim1];
#line 599 "ssytri_3x.f"
		    }
#line 600 "ssytri_3x.f"
		}

#line 602 "ssytri_3x.f"
	    } else {

/*              L11 =  L11**T * invD1 * L11 */

#line 606 "ssytri_3x.f"
		i__1 = nnb;
#line 606 "ssytri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 607 "ssytri_3x.f"
		    i__2 = i__;
#line 607 "ssytri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 608 "ssytri_3x.f"
			a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + 
				j * work_dim1];
#line 609 "ssytri_3x.f"
		    }
#line 610 "ssytri_3x.f"
		}
#line 611 "ssytri_3x.f"
	    }

/*           Next Block */

#line 615 "ssytri_3x.f"
	    cut += nnb;

#line 617 "ssytri_3x.f"
	}

/*        Apply PERMUTATIONS P and P**T: */
/*        P * inv(L**T) * inv(D) * inv(L) * P**T. */
/*        Interchange rows and columns I and IPIV(I) in reverse order */
/*        from the formation order of IPIV vector for Lower case. */

/*        ( We can use a loop over IPIV with increment -1, */
/*        since the ABS value of IPIV(I) represents the row (column) */
/*        index of the interchange with row (column) i in both 1x1 */
/*        and 2x2 pivot cases, i.e. we don't need separate code branches */
/*        for 1x1 and 2x2 pivot cases ) */

#line 630 "ssytri_3x.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 631 "ssytri_3x.f"
	    ip = (i__1 = ipiv[i__], abs(i__1));
#line 632 "ssytri_3x.f"
	    if (ip != i__) {
#line 633 "ssytri_3x.f"
		if (i__ < ip) {
#line 633 "ssytri_3x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 633 "ssytri_3x.f"
		}
#line 634 "ssytri_3x.f"
		if (i__ > ip) {
#line 634 "ssytri_3x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 634 "ssytri_3x.f"
		}
#line 635 "ssytri_3x.f"
	    }
#line 636 "ssytri_3x.f"
	}

#line 638 "ssytri_3x.f"
    }

#line 640 "ssytri_3x.f"
    return 0;

/*     End of SSYTRI_3X */

} /* ssytri_3x__ */

