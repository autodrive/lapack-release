#line 1 "chetri_3x.f"
/* chetri_3x.f -- translated by f2c (version 20100827).
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

#line 1 "chetri_3x.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};

/* > \brief \b CHETRI_3X */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRI_3X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri_
3x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri_
3x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri_
3x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ),  E( * ), WORK( N+NB+1, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > CHETRI_3X computes the inverse of a complex Hermitian indefinite */
/* > matrix A using the factorization computed by CHETRF_RK or CHETRF_BK: */
/* > */
/* >     A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is Hermitian and block */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, diagonal of the block diagonal matrix D and */
/* >          factors U or L as computed by CHETRF_RK and CHETRF_BK: */
/* >            a) ONLY diagonal elements of the Hermitian block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                should be provided on entry in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, if INFO = 0, the Hermitian inverse of the original */
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
/* >          E is COMPLEX array, dimension (N) */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the Hermitian block diagonal matrix D */
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
/* >          as determined by CHETRF_RK or CHETRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N+NB+1,NB+3). */
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

/* > \ingroup complexHEcomputational */

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
/* Subroutine */ int chetri_3x__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *e, integer *ipiv, doublecomplex *work, 
	integer *nb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    extern /* Subroutine */ int cheswapr_(char *, integer *, doublecomplex *, 
	    integer *, integer *, integer *, ftnlen);
    static doublecomplex d__;
    static integer i__, j, k;
    static doublereal t, ak;
    static integer u11, ip, nnb, cut;
    static doublereal akp1;
    static integer invd;
    static doublecomplex akkp1;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    static doublecomplex u01_i_j__, u11_i_j__;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer icount;
    extern /* Subroutine */ int ctrtri_(char *, char *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);
    static doublecomplex u01_ip1_j__, u11_ip1_j__;


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

#line 206 "chetri_3x.f"
    /* Parameter adjustments */
#line 206 "chetri_3x.f"
    a_dim1 = *lda;
#line 206 "chetri_3x.f"
    a_offset = 1 + a_dim1;
#line 206 "chetri_3x.f"
    a -= a_offset;
#line 206 "chetri_3x.f"
    --e;
#line 206 "chetri_3x.f"
    --ipiv;
#line 206 "chetri_3x.f"
    work_dim1 = *n + *nb + 1;
#line 206 "chetri_3x.f"
    work_offset = 1 + work_dim1;
#line 206 "chetri_3x.f"
    work -= work_offset;
#line 206 "chetri_3x.f"

#line 206 "chetri_3x.f"
    /* Function Body */
#line 206 "chetri_3x.f"
    *info = 0;
#line 207 "chetri_3x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 208 "chetri_3x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 209 "chetri_3x.f"
	*info = -1;
#line 210 "chetri_3x.f"
    } else if (*n < 0) {
#line 211 "chetri_3x.f"
	*info = -2;
#line 212 "chetri_3x.f"
    } else if (*lda < max(1,*n)) {
#line 213 "chetri_3x.f"
	*info = -4;
#line 214 "chetri_3x.f"
    }

/*     Quick return if possible */

#line 218 "chetri_3x.f"
    if (*info != 0) {
#line 219 "chetri_3x.f"
	i__1 = -(*info);
#line 219 "chetri_3x.f"
	xerbla_("CHETRI_3X", &i__1, (ftnlen)9);
#line 220 "chetri_3x.f"
	return 0;
#line 221 "chetri_3x.f"
    }
#line 222 "chetri_3x.f"
    if (*n == 0) {
#line 222 "chetri_3x.f"
	return 0;
#line 222 "chetri_3x.f"
    }

/*     Workspace got Non-diag elements of D */

#line 227 "chetri_3x.f"
    i__1 = *n;
#line 227 "chetri_3x.f"
    for (k = 1; k <= i__1; ++k) {
#line 228 "chetri_3x.f"
	i__2 = k + work_dim1;
#line 228 "chetri_3x.f"
	i__3 = k;
#line 228 "chetri_3x.f"
	work[i__2].r = e[i__3].r, work[i__2].i = e[i__3].i;
#line 229 "chetri_3x.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 233 "chetri_3x.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 237 "chetri_3x.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 238 "chetri_3x.f"
	    i__1 = *info + *info * a_dim1;
#line 238 "chetri_3x.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 238 "chetri_3x.f"
		return 0;
#line 238 "chetri_3x.f"
	    }
#line 240 "chetri_3x.f"
	}
#line 241 "chetri_3x.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 245 "chetri_3x.f"
	i__1 = *n;
#line 245 "chetri_3x.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 246 "chetri_3x.f"
	    i__2 = *info + *info * a_dim1;
#line 246 "chetri_3x.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 246 "chetri_3x.f"
		return 0;
#line 246 "chetri_3x.f"
	    }
#line 248 "chetri_3x.f"
	}
#line 249 "chetri_3x.f"
    }

#line 251 "chetri_3x.f"
    *info = 0;

/*     Splitting Workspace */
/*     U01 is a block ( N, NB+1 ) */
/*     The first element of U01 is in WORK( 1, 1 ) */
/*     U11 is a block ( NB+1, NB+1 ) */
/*     The first element of U11 is in WORK( N+1, 1 ) */

#line 259 "chetri_3x.f"
    u11 = *n;

/*     INVD is a block ( N, 2 ) */
/*     The first element of INVD is in WORK( 1, INVD ) */

#line 264 "chetri_3x.f"
    invd = *nb + 2;
#line 266 "chetri_3x.f"
    if (upper) {

/*        Begin Upper */

/*        invA = P * inv(U**H) * inv(D) * inv(U) * P**T. */

#line 272 "chetri_3x.f"
	ctrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*        inv(D) and inv(D) * inv(U) */

#line 276 "chetri_3x.f"
	k = 1;
#line 277 "chetri_3x.f"
	while(k <= *n) {
#line 278 "chetri_3x.f"
	    if (ipiv[k] > 0) {
/*              1 x 1 diagonal NNB */
#line 280 "chetri_3x.f"
		i__1 = k + invd * work_dim1;
#line 280 "chetri_3x.f"
		i__2 = k + k * a_dim1;
#line 280 "chetri_3x.f"
		d__1 = 1. / a[i__2].r;
#line 280 "chetri_3x.f"
		work[i__1].r = d__1, work[i__1].i = 0.;
#line 281 "chetri_3x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 281 "chetri_3x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 282 "chetri_3x.f"
	    } else {
/*              2 x 2 diagonal NNB */
#line 284 "chetri_3x.f"
		t = z_abs(&work[k + 1 + work_dim1]);
#line 285 "chetri_3x.f"
		i__1 = k + k * a_dim1;
#line 285 "chetri_3x.f"
		ak = a[i__1].r / t;
#line 286 "chetri_3x.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 286 "chetri_3x.f"
		akp1 = a[i__1].r / t;
#line 287 "chetri_3x.f"
		i__1 = k + 1 + work_dim1;
#line 287 "chetri_3x.f"
		z__1.r = work[i__1].r / t, z__1.i = work[i__1].i / t;
#line 287 "chetri_3x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 288 "chetri_3x.f"
		d__1 = ak * akp1;
#line 288 "chetri_3x.f"
		z__2.r = d__1 - 1., z__2.i = -0.;
#line 288 "chetri_3x.f"
		z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 288 "chetri_3x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 289 "chetri_3x.f"
		i__1 = k + invd * work_dim1;
#line 289 "chetri_3x.f"
		z__2.r = akp1, z__2.i = 0.;
#line 289 "chetri_3x.f"
		z_div(&z__1, &z__2, &d__);
#line 289 "chetri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 290 "chetri_3x.f"
		i__1 = k + 1 + (invd + 1) * work_dim1;
#line 290 "chetri_3x.f"
		z__2.r = ak, z__2.i = 0.;
#line 290 "chetri_3x.f"
		z_div(&z__1, &z__2, &d__);
#line 290 "chetri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 291 "chetri_3x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 291 "chetri_3x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 291 "chetri_3x.f"
		z_div(&z__1, &z__2, &d__);
#line 291 "chetri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 292 "chetri_3x.f"
		i__1 = k + 1 + invd * work_dim1;
#line 292 "chetri_3x.f"
		d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
#line 292 "chetri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 293 "chetri_3x.f"
		++k;
#line 294 "chetri_3x.f"
	    }
#line 295 "chetri_3x.f"
	    ++k;
#line 296 "chetri_3x.f"
	}

/*        inv(U**H) = (inv(U))**H */

/*        inv(U**H) * inv(D) * inv(U) */

#line 302 "chetri_3x.f"
	cut = *n;
#line 303 "chetri_3x.f"
	while(cut > 0) {
#line 304 "chetri_3x.f"
	    nnb = *nb;
#line 305 "chetri_3x.f"
	    if (cut <= nnb) {
#line 306 "chetri_3x.f"
		nnb = cut;
#line 307 "chetri_3x.f"
	    } else {
#line 308 "chetri_3x.f"
		icount = 0;
/*              count negative elements, */
#line 310 "chetri_3x.f"
		i__1 = cut;
#line 310 "chetri_3x.f"
		for (i__ = cut + 1 - nnb; i__ <= i__1; ++i__) {
#line 311 "chetri_3x.f"
		    if (ipiv[i__] < 0) {
#line 311 "chetri_3x.f"
			++icount;
#line 311 "chetri_3x.f"
		    }
#line 312 "chetri_3x.f"
		}
/*              need a even number for a clear cut */
#line 314 "chetri_3x.f"
		if (icount % 2 == 1) {
#line 314 "chetri_3x.f"
		    ++nnb;
#line 314 "chetri_3x.f"
		}
#line 315 "chetri_3x.f"
	    }
#line 317 "chetri_3x.f"
	    cut -= nnb;

/*           U01 Block */

#line 321 "chetri_3x.f"
	    i__1 = cut;
#line 321 "chetri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 322 "chetri_3x.f"
		i__2 = nnb;
#line 322 "chetri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 323 "chetri_3x.f"
		    i__3 = i__ + j * work_dim1;
#line 323 "chetri_3x.f"
		    i__4 = i__ + (cut + j) * a_dim1;
#line 323 "chetri_3x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 324 "chetri_3x.f"
		}
#line 325 "chetri_3x.f"
	    }

/*           U11 Block */

#line 329 "chetri_3x.f"
	    i__1 = nnb;
#line 329 "chetri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 330 "chetri_3x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 330 "chetri_3x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 331 "chetri_3x.f"
		i__2 = i__ - 1;
#line 331 "chetri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 332 "chetri_3x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 332 "chetri_3x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 333 "chetri_3x.f"
		}
#line 334 "chetri_3x.f"
		i__2 = nnb;
#line 334 "chetri_3x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 335 "chetri_3x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 335 "chetri_3x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 335 "chetri_3x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 336 "chetri_3x.f"
		}
#line 337 "chetri_3x.f"
	    }

/*           invD * U01 */

#line 341 "chetri_3x.f"
	    i__ = 1;
#line 342 "chetri_3x.f"
	    while(i__ <= cut) {
#line 343 "chetri_3x.f"
		if (ipiv[i__] > 0) {
#line 344 "chetri_3x.f"
		    i__1 = nnb;
#line 344 "chetri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 345 "chetri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 345 "chetri_3x.f"
			i__3 = i__ + invd * work_dim1;
#line 345 "chetri_3x.f"
			i__4 = i__ + j * work_dim1;
#line 345 "chetri_3x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 345 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 346 "chetri_3x.f"
		    }
#line 347 "chetri_3x.f"
		} else {
#line 348 "chetri_3x.f"
		    i__1 = nnb;
#line 348 "chetri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 349 "chetri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 349 "chetri_3x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 350 "chetri_3x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 350 "chetri_3x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 351 "chetri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 351 "chetri_3x.f"
			i__3 = i__ + invd * work_dim1;
#line 351 "chetri_3x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 351 "chetri_3x.f"
			i__4 = i__ + (invd + 1) * work_dim1;
#line 351 "chetri_3x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 351 "chetri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 351 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 353 "chetri_3x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 353 "chetri_3x.f"
			i__3 = i__ + 1 + invd * work_dim1;
#line 353 "chetri_3x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 353 "chetri_3x.f"
			i__4 = i__ + 1 + (invd + 1) * work_dim1;
#line 353 "chetri_3x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 353 "chetri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 353 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 355 "chetri_3x.f"
		    }
#line 356 "chetri_3x.f"
		    ++i__;
#line 357 "chetri_3x.f"
		}
#line 358 "chetri_3x.f"
		++i__;
#line 359 "chetri_3x.f"
	    }

/*           invD1 * U11 */

#line 363 "chetri_3x.f"
	    i__ = 1;
#line 364 "chetri_3x.f"
	    while(i__ <= nnb) {
#line 365 "chetri_3x.f"
		if (ipiv[cut + i__] > 0) {
#line 366 "chetri_3x.f"
		    i__1 = nnb;
#line 366 "chetri_3x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 367 "chetri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 367 "chetri_3x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 367 "chetri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 367 "chetri_3x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 367 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 368 "chetri_3x.f"
		    }
#line 369 "chetri_3x.f"
		} else {
#line 370 "chetri_3x.f"
		    i__1 = nnb;
#line 370 "chetri_3x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 371 "chetri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 371 "chetri_3x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 372 "chetri_3x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 372 "chetri_3x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 373 "chetri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 373 "chetri_3x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 373 "chetri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 373 "chetri_3x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 373 "chetri_3x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 373 "chetri_3x.f"
			i__6 = u11 + i__ + 1 + j * work_dim1;
#line 373 "chetri_3x.f"
			z__3.r = work[i__5].r * work[i__6].r - work[i__5].i * 
				work[i__6].i, z__3.i = work[i__5].r * work[
				i__6].i + work[i__5].i * work[i__6].r;
#line 373 "chetri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 373 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 375 "chetri_3x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 375 "chetri_3x.f"
			i__3 = cut + i__ + 1 + invd * work_dim1;
#line 375 "chetri_3x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 375 "chetri_3x.f"
			i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
#line 375 "chetri_3x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 375 "chetri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 375 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 377 "chetri_3x.f"
		    }
#line 378 "chetri_3x.f"
		    ++i__;
#line 379 "chetri_3x.f"
		}
#line 380 "chetri_3x.f"
		++i__;
#line 381 "chetri_3x.f"
	    }

/*           U11**H * invD1 * U11 -> U11 */

#line 385 "chetri_3x.f"
	    i__1 = *n + *nb + 1;
#line 385 "chetri_3x.f"
	    ctrmm_("L", "U", "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 
		    1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 389 "chetri_3x.f"
	    i__1 = nnb;
#line 389 "chetri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 390 "chetri_3x.f"
		i__2 = nnb;
#line 390 "chetri_3x.f"
		for (j = i__; j <= i__2; ++j) {
#line 391 "chetri_3x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 391 "chetri_3x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 391 "chetri_3x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 392 "chetri_3x.f"
		}
#line 393 "chetri_3x.f"
	    }

/*           U01**H * invD * U01 -> A( CUT+I, CUT+J ) */

#line 397 "chetri_3x.f"
	    i__1 = *n + *nb + 1;
#line 397 "chetri_3x.f"
	    i__2 = *n + *nb + 1;
#line 397 "chetri_3x.f"
	    cgemm_("C", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 
		    1], lda, &work[work_offset], &i__1, &c_b2, &work[u11 + 1 
		    + work_dim1], &i__2, (ftnlen)1, (ftnlen)1);

/*           U11 =  U11**H * invD1 * U11 + U01**H * invD * U01 */

#line 404 "chetri_3x.f"
	    i__1 = nnb;
#line 404 "chetri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 405 "chetri_3x.f"
		i__2 = nnb;
#line 405 "chetri_3x.f"
		for (j = i__; j <= i__2; ++j) {
#line 406 "chetri_3x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 406 "chetri_3x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 406 "chetri_3x.f"
		    i__5 = u11 + i__ + j * work_dim1;
#line 406 "chetri_3x.f"
		    z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i + 
			    work[i__5].i;
#line 406 "chetri_3x.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 407 "chetri_3x.f"
		}
#line 408 "chetri_3x.f"
	    }

/*           U01 =  U00**H * invD0 * U01 */

#line 412 "chetri_3x.f"
	    i__1 = *n + *nb + 1;
#line 412 "chetri_3x.f"
	    ctrmm_("L", uplo, "C", "U", &cut, &nnb, &c_b1, &a[a_offset], lda, 
		    &work[work_offset], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);

/*           Update U01 */

#line 418 "chetri_3x.f"
	    i__1 = cut;
#line 418 "chetri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 419 "chetri_3x.f"
		i__2 = nnb;
#line 419 "chetri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 420 "chetri_3x.f"
		    i__3 = i__ + (cut + j) * a_dim1;
#line 420 "chetri_3x.f"
		    i__4 = i__ + j * work_dim1;
#line 420 "chetri_3x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 421 "chetri_3x.f"
		}
#line 422 "chetri_3x.f"
	    }

/*           Next Block */

#line 426 "chetri_3x.f"
	}

/*        Apply PERMUTATIONS P and P**T: */
/*        P * inv(U**H) * inv(D) * inv(U) * P**T. */
/*        Interchange rows and columns I and IPIV(I) in reverse order */
/*        from the formation order of IPIV vector for Upper case. */

/*        ( We can use a loop over IPIV with increment 1, */
/*        since the ABS value of IPIV(I) represents the row (column) */
/*        index of the interchange with row (column) i in both 1x1 */
/*        and 2x2 pivot cases, i.e. we don't need separate code branches */
/*        for 1x1 and 2x2 pivot cases ) */

#line 439 "chetri_3x.f"
	i__1 = *n;
#line 439 "chetri_3x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 440 "chetri_3x.f"
	    ip = (i__2 = ipiv[i__], abs(i__2));
#line 441 "chetri_3x.f"
	    if (ip != i__) {
#line 442 "chetri_3x.f"
		if (i__ < ip) {
#line 442 "chetri_3x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 442 "chetri_3x.f"
		}
#line 443 "chetri_3x.f"
		if (i__ > ip) {
#line 443 "chetri_3x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 443 "chetri_3x.f"
		}
#line 444 "chetri_3x.f"
	    }
#line 445 "chetri_3x.f"
	}

#line 447 "chetri_3x.f"
    } else {

/*        Begin Lower */

/*        inv A = P * inv(L**H) * inv(D) * inv(L) * P**T. */

#line 453 "chetri_3x.f"
	ctrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*        inv(D) and inv(D) * inv(L) */

#line 457 "chetri_3x.f"
	k = *n;
#line 458 "chetri_3x.f"
	while(k >= 1) {
#line 459 "chetri_3x.f"
	    if (ipiv[k] > 0) {
/*              1 x 1 diagonal NNB */
#line 461 "chetri_3x.f"
		i__1 = k + invd * work_dim1;
#line 461 "chetri_3x.f"
		i__2 = k + k * a_dim1;
#line 461 "chetri_3x.f"
		d__1 = 1. / a[i__2].r;
#line 461 "chetri_3x.f"
		work[i__1].r = d__1, work[i__1].i = 0.;
#line 462 "chetri_3x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 462 "chetri_3x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 463 "chetri_3x.f"
	    } else {
/*              2 x 2 diagonal NNB */
#line 465 "chetri_3x.f"
		t = z_abs(&work[k - 1 + work_dim1]);
#line 466 "chetri_3x.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 466 "chetri_3x.f"
		ak = a[i__1].r / t;
#line 467 "chetri_3x.f"
		i__1 = k + k * a_dim1;
#line 467 "chetri_3x.f"
		akp1 = a[i__1].r / t;
#line 468 "chetri_3x.f"
		i__1 = k - 1 + work_dim1;
#line 468 "chetri_3x.f"
		z__1.r = work[i__1].r / t, z__1.i = work[i__1].i / t;
#line 468 "chetri_3x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 469 "chetri_3x.f"
		d__1 = ak * akp1;
#line 469 "chetri_3x.f"
		z__2.r = d__1 - 1., z__2.i = -0.;
#line 469 "chetri_3x.f"
		z__1.r = t * z__2.r, z__1.i = t * z__2.i;
#line 469 "chetri_3x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 470 "chetri_3x.f"
		i__1 = k - 1 + invd * work_dim1;
#line 470 "chetri_3x.f"
		z__2.r = akp1, z__2.i = 0.;
#line 470 "chetri_3x.f"
		z_div(&z__1, &z__2, &d__);
#line 470 "chetri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 471 "chetri_3x.f"
		i__1 = k + invd * work_dim1;
#line 471 "chetri_3x.f"
		z__2.r = ak, z__2.i = 0.;
#line 471 "chetri_3x.f"
		z_div(&z__1, &z__2, &d__);
#line 471 "chetri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 472 "chetri_3x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 472 "chetri_3x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 472 "chetri_3x.f"
		z_div(&z__1, &z__2, &d__);
#line 472 "chetri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 473 "chetri_3x.f"
		i__1 = k - 1 + (invd + 1) * work_dim1;
#line 473 "chetri_3x.f"
		d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
#line 473 "chetri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 474 "chetri_3x.f"
		--k;
#line 475 "chetri_3x.f"
	    }
#line 476 "chetri_3x.f"
	    --k;
#line 477 "chetri_3x.f"
	}

/*        inv(L**H) = (inv(L))**H */

/*        inv(L**H) * inv(D) * inv(L) */

#line 483 "chetri_3x.f"
	cut = 0;
#line 484 "chetri_3x.f"
	while(cut < *n) {
#line 485 "chetri_3x.f"
	    nnb = *nb;
#line 486 "chetri_3x.f"
	    if (cut + nnb > *n) {
#line 487 "chetri_3x.f"
		nnb = *n - cut;
#line 488 "chetri_3x.f"
	    } else {
#line 489 "chetri_3x.f"
		icount = 0;
/*              count negative elements, */
#line 491 "chetri_3x.f"
		i__1 = cut + nnb;
#line 491 "chetri_3x.f"
		for (i__ = cut + 1; i__ <= i__1; ++i__) {
#line 492 "chetri_3x.f"
		    if (ipiv[i__] < 0) {
#line 492 "chetri_3x.f"
			++icount;
#line 492 "chetri_3x.f"
		    }
#line 493 "chetri_3x.f"
		}
/*              need a even number for a clear cut */
#line 495 "chetri_3x.f"
		if (icount % 2 == 1) {
#line 495 "chetri_3x.f"
		    ++nnb;
#line 495 "chetri_3x.f"
		}
#line 496 "chetri_3x.f"
	    }

/*           L21 Block */

#line 500 "chetri_3x.f"
	    i__1 = *n - cut - nnb;
#line 500 "chetri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 501 "chetri_3x.f"
		i__2 = nnb;
#line 501 "chetri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 502 "chetri_3x.f"
		    i__3 = i__ + j * work_dim1;
#line 502 "chetri_3x.f"
		    i__4 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 502 "chetri_3x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 503 "chetri_3x.f"
		}
#line 504 "chetri_3x.f"
	    }

/*           L11 Block */

#line 508 "chetri_3x.f"
	    i__1 = nnb;
#line 508 "chetri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 509 "chetri_3x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 509 "chetri_3x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 510 "chetri_3x.f"
		i__2 = nnb;
#line 510 "chetri_3x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 511 "chetri_3x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 511 "chetri_3x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 512 "chetri_3x.f"
		}
#line 513 "chetri_3x.f"
		i__2 = i__ - 1;
#line 513 "chetri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 514 "chetri_3x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 514 "chetri_3x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 514 "chetri_3x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 515 "chetri_3x.f"
		}
#line 516 "chetri_3x.f"
	    }

/*           invD*L21 */

#line 520 "chetri_3x.f"
	    i__ = *n - cut - nnb;
#line 521 "chetri_3x.f"
	    while(i__ >= 1) {
#line 522 "chetri_3x.f"
		if (ipiv[cut + nnb + i__] > 0) {
#line 523 "chetri_3x.f"
		    i__1 = nnb;
#line 523 "chetri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 524 "chetri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 524 "chetri_3x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 524 "chetri_3x.f"
			i__4 = i__ + j * work_dim1;
#line 524 "chetri_3x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 524 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 525 "chetri_3x.f"
		    }
#line 526 "chetri_3x.f"
		} else {
#line 527 "chetri_3x.f"
		    i__1 = nnb;
#line 527 "chetri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 528 "chetri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 528 "chetri_3x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 529 "chetri_3x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 529 "chetri_3x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 530 "chetri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 530 "chetri_3x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 530 "chetri_3x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 530 "chetri_3x.f"
			i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
#line 530 "chetri_3x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 530 "chetri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 530 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 532 "chetri_3x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 532 "chetri_3x.f"
			i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
#line 532 "chetri_3x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 532 "chetri_3x.f"
			i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
#line 532 "chetri_3x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 532 "chetri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 532 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 534 "chetri_3x.f"
		    }
#line 535 "chetri_3x.f"
		    --i__;
#line 536 "chetri_3x.f"
		}
#line 537 "chetri_3x.f"
		--i__;
#line 538 "chetri_3x.f"
	    }

/*           invD1*L11 */

#line 542 "chetri_3x.f"
	    i__ = nnb;
#line 543 "chetri_3x.f"
	    while(i__ >= 1) {
#line 544 "chetri_3x.f"
		if (ipiv[cut + i__] > 0) {
#line 545 "chetri_3x.f"
		    i__1 = nnb;
#line 545 "chetri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 546 "chetri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 546 "chetri_3x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 546 "chetri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 546 "chetri_3x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 546 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 547 "chetri_3x.f"
		    }
#line 549 "chetri_3x.f"
		} else {
#line 550 "chetri_3x.f"
		    i__1 = nnb;
#line 550 "chetri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 551 "chetri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 551 "chetri_3x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 552 "chetri_3x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 552 "chetri_3x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 553 "chetri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 553 "chetri_3x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 553 "chetri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 553 "chetri_3x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 553 "chetri_3x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 553 "chetri_3x.f"
			z__3.r = work[i__5].r * u11_ip1_j__.r - work[i__5].i *
				 u11_ip1_j__.i, z__3.i = work[i__5].r * 
				u11_ip1_j__.i + work[i__5].i * u11_ip1_j__.r;
#line 553 "chetri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 553 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 555 "chetri_3x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 555 "chetri_3x.f"
			i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
#line 555 "chetri_3x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 555 "chetri_3x.f"
			i__4 = cut + i__ - 1 + invd * work_dim1;
#line 555 "chetri_3x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 555 "chetri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 555 "chetri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 557 "chetri_3x.f"
		    }
#line 558 "chetri_3x.f"
		    --i__;
#line 559 "chetri_3x.f"
		}
#line 560 "chetri_3x.f"
		--i__;
#line 561 "chetri_3x.f"
	    }

/*           L11**H * invD1 * L11 -> L11 */

#line 565 "chetri_3x.f"
	    i__1 = *n + *nb + 1;
#line 565 "chetri_3x.f"
	    ctrmm_("L", uplo, "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 570 "chetri_3x.f"
	    i__1 = nnb;
#line 570 "chetri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 571 "chetri_3x.f"
		i__2 = i__;
#line 571 "chetri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 572 "chetri_3x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 572 "chetri_3x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 572 "chetri_3x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 573 "chetri_3x.f"
		}
#line 574 "chetri_3x.f"
	    }

#line 576 "chetri_3x.f"
	    if (cut + nnb < *n) {

/*              L21**H * invD2*L21 -> A( CUT+I, CUT+J ) */

#line 580 "chetri_3x.f"
		i__1 = *n - nnb - cut;
#line 580 "chetri_3x.f"
		i__2 = *n + *nb + 1;
#line 580 "chetri_3x.f"
		i__3 = *n + *nb + 1;
#line 580 "chetri_3x.f"
		cgemm_("C", "N", &nnb, &nnb, &i__1, &c_b1, &a[cut + nnb + 1 + 
			(cut + 1) * a_dim1], lda, &work[work_offset], &i__2, &
			c_b2, &work[u11 + 1 + work_dim1], &i__3, (ftnlen)1, (
			ftnlen)1);

/*              L11 =  L11**H * invD1 * L11 + U01**H * invD * U01 */

#line 587 "chetri_3x.f"
		i__1 = nnb;
#line 587 "chetri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 588 "chetri_3x.f"
		    i__2 = i__;
#line 588 "chetri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 589 "chetri_3x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 589 "chetri_3x.f"
			i__4 = cut + i__ + (cut + j) * a_dim1;
#line 589 "chetri_3x.f"
			i__5 = u11 + i__ + j * work_dim1;
#line 589 "chetri_3x.f"
			z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i 
				+ work[i__5].i;
#line 589 "chetri_3x.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 590 "chetri_3x.f"
		    }
#line 591 "chetri_3x.f"
		}

/*              L01 =  L22**H * invD2 * L21 */

#line 595 "chetri_3x.f"
		i__1 = *n - nnb - cut;
#line 595 "chetri_3x.f"
		i__2 = *n + *nb + 1;
#line 595 "chetri_3x.f"
		ctrmm_("L", uplo, "C", "U", &i__1, &nnb, &c_b1, &a[cut + nnb 
			+ 1 + (cut + nnb + 1) * a_dim1], lda, &work[
			work_offset], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);

/*              Update L21 */

#line 601 "chetri_3x.f"
		i__1 = *n - cut - nnb;
#line 601 "chetri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 602 "chetri_3x.f"
		    i__2 = nnb;
#line 602 "chetri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 603 "chetri_3x.f"
			i__3 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 603 "chetri_3x.f"
			i__4 = i__ + j * work_dim1;
#line 603 "chetri_3x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 604 "chetri_3x.f"
		    }
#line 605 "chetri_3x.f"
		}

#line 607 "chetri_3x.f"
	    } else {

/*              L11 =  L11**H * invD1 * L11 */

#line 611 "chetri_3x.f"
		i__1 = nnb;
#line 611 "chetri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 612 "chetri_3x.f"
		    i__2 = i__;
#line 612 "chetri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 613 "chetri_3x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 613 "chetri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 613 "chetri_3x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 614 "chetri_3x.f"
		    }
#line 615 "chetri_3x.f"
		}
#line 616 "chetri_3x.f"
	    }

/*           Next Block */

#line 620 "chetri_3x.f"
	    cut += nnb;

#line 622 "chetri_3x.f"
	}

/*        Apply PERMUTATIONS P and P**T: */
/*        P * inv(L**H) * inv(D) * inv(L) * P**T. */
/*        Interchange rows and columns I and IPIV(I) in reverse order */
/*        from the formation order of IPIV vector for Lower case. */

/*        ( We can use a loop over IPIV with increment -1, */
/*        since the ABS value of IPIV(I) represents the row (column) */
/*        index of the interchange with row (column) i in both 1x1 */
/*        and 2x2 pivot cases, i.e. we don't need separate code branches */
/*        for 1x1 and 2x2 pivot cases ) */

#line 635 "chetri_3x.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 636 "chetri_3x.f"
	    ip = (i__1 = ipiv[i__], abs(i__1));
#line 637 "chetri_3x.f"
	    if (ip != i__) {
#line 638 "chetri_3x.f"
		if (i__ < ip) {
#line 638 "chetri_3x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 638 "chetri_3x.f"
		}
#line 639 "chetri_3x.f"
		if (i__ > ip) {
#line 639 "chetri_3x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 639 "chetri_3x.f"
		}
#line 640 "chetri_3x.f"
	    }
#line 641 "chetri_3x.f"
	}

#line 643 "chetri_3x.f"
    }

#line 645 "chetri_3x.f"
    return 0;

/*     End of CHETRI_3X */

} /* chetri_3x__ */

