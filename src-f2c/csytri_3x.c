#line 1 "csytri_3x.f"
/* csytri_3x.f -- translated by f2c (version 20100827).
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

#line 1 "csytri_3x.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};

/* > \brief \b CSYTRI_3X */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRI_3X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytri_
3x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytri_
3x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytri_
3x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ),  E( * ), WORK( N+NB+1, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > CSYTRI_3X computes the inverse of a complex symmetric indefinite */
/* > matrix A using the factorization computed by CSYTRF_RK or CSYTRF_BK: */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, diagonal of the block diagonal matrix D and */
/* >          factors U or L as computed by CSYTRF_RK and CSYTRF_BK: */
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
/* >          E is COMPLEX array, dimension (N) */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) not refernced; */
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
/* >          as determined by CSYTRF_RK or CSYTRF_BK. */
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

/* > \date December 2016 */

/* > \ingroup complexSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > \verbatim */
/* > */
/* >  December 2016,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int csytri_3x__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *e, integer *ipiv, doublecomplex *work, 
	integer *nb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    extern /* Subroutine */ int csyswapr_(char *, integer *, doublecomplex *, 
	    integer *, integer *, integer *, ftnlen);
    static doublecomplex d__;
    static integer i__, j, k;
    static doublecomplex t, ak;
    static integer u11, ip, nnb, cut;
    static doublecomplex akp1;
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

#line 203 "csytri_3x.f"
    /* Parameter adjustments */
#line 203 "csytri_3x.f"
    a_dim1 = *lda;
#line 203 "csytri_3x.f"
    a_offset = 1 + a_dim1;
#line 203 "csytri_3x.f"
    a -= a_offset;
#line 203 "csytri_3x.f"
    --e;
#line 203 "csytri_3x.f"
    --ipiv;
#line 203 "csytri_3x.f"
    work_dim1 = *n + *nb + 1;
#line 203 "csytri_3x.f"
    work_offset = 1 + work_dim1;
#line 203 "csytri_3x.f"
    work -= work_offset;
#line 203 "csytri_3x.f"

#line 203 "csytri_3x.f"
    /* Function Body */
#line 203 "csytri_3x.f"
    *info = 0;
#line 204 "csytri_3x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 205 "csytri_3x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 206 "csytri_3x.f"
	*info = -1;
#line 207 "csytri_3x.f"
    } else if (*n < 0) {
#line 208 "csytri_3x.f"
	*info = -2;
#line 209 "csytri_3x.f"
    } else if (*lda < max(1,*n)) {
#line 210 "csytri_3x.f"
	*info = -4;
#line 211 "csytri_3x.f"
    }

/*     Quick return if possible */

#line 215 "csytri_3x.f"
    if (*info != 0) {
#line 216 "csytri_3x.f"
	i__1 = -(*info);
#line 216 "csytri_3x.f"
	xerbla_("CSYTRI_3X", &i__1, (ftnlen)9);
#line 217 "csytri_3x.f"
	return 0;
#line 218 "csytri_3x.f"
    }
#line 219 "csytri_3x.f"
    if (*n == 0) {
#line 219 "csytri_3x.f"
	return 0;
#line 219 "csytri_3x.f"
    }

/*     Workspace got Non-diag elements of D */

#line 224 "csytri_3x.f"
    i__1 = *n;
#line 224 "csytri_3x.f"
    for (k = 1; k <= i__1; ++k) {
#line 225 "csytri_3x.f"
	i__2 = k + work_dim1;
#line 225 "csytri_3x.f"
	i__3 = k;
#line 225 "csytri_3x.f"
	work[i__2].r = e[i__3].r, work[i__2].i = e[i__3].i;
#line 226 "csytri_3x.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 230 "csytri_3x.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 234 "csytri_3x.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 235 "csytri_3x.f"
	    i__1 = *info + *info * a_dim1;
#line 235 "csytri_3x.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 235 "csytri_3x.f"
		return 0;
#line 235 "csytri_3x.f"
	    }
#line 237 "csytri_3x.f"
	}
#line 238 "csytri_3x.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 242 "csytri_3x.f"
	i__1 = *n;
#line 242 "csytri_3x.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 243 "csytri_3x.f"
	    i__2 = *info + *info * a_dim1;
#line 243 "csytri_3x.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 243 "csytri_3x.f"
		return 0;
#line 243 "csytri_3x.f"
	    }
#line 245 "csytri_3x.f"
	}
#line 246 "csytri_3x.f"
    }

#line 248 "csytri_3x.f"
    *info = 0;

/*     Splitting Workspace */
/*     U01 is a block ( N, NB+1 ) */
/*     The first element of U01 is in WORK( 1, 1 ) */
/*     U11 is a block ( NB+1, NB+1 ) */
/*     The first element of U11 is in WORK( N+1, 1 ) */

#line 256 "csytri_3x.f"
    u11 = *n;

/*     INVD is a block ( N, 2 ) */
/*     The first element of INVD is in WORK( 1, INVD ) */

#line 261 "csytri_3x.f"
    invd = *nb + 2;
#line 263 "csytri_3x.f"
    if (upper) {

/*        Begin Upper */

/*        invA = P * inv(U**T) * inv(D) * inv(U) * P**T. */

#line 269 "csytri_3x.f"
	ctrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*        inv(D) and inv(D) * inv(U) */

#line 273 "csytri_3x.f"
	k = 1;
#line 274 "csytri_3x.f"
	while(k <= *n) {
#line 275 "csytri_3x.f"
	    if (ipiv[k] > 0) {
/*              1 x 1 diagonal NNB */
#line 277 "csytri_3x.f"
		i__1 = k + invd * work_dim1;
#line 277 "csytri_3x.f"
		z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 277 "csytri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 278 "csytri_3x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 278 "csytri_3x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 279 "csytri_3x.f"
	    } else {
/*              2 x 2 diagonal NNB */
#line 281 "csytri_3x.f"
		i__1 = k + 1 + work_dim1;
#line 281 "csytri_3x.f"
		t.r = work[i__1].r, t.i = work[i__1].i;
#line 282 "csytri_3x.f"
		z_div(&z__1, &a[k + k * a_dim1], &t);
#line 282 "csytri_3x.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 283 "csytri_3x.f"
		z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &t);
#line 283 "csytri_3x.f"
		akp1.r = z__1.r, akp1.i = z__1.i;
#line 284 "csytri_3x.f"
		z_div(&z__1, &work[k + 1 + work_dim1], &t);
#line 284 "csytri_3x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 285 "csytri_3x.f"
		z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * 
			akp1.i + ak.i * akp1.r;
#line 285 "csytri_3x.f"
		z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 285 "csytri_3x.f"
		z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + 
			t.i * z__2.r;
#line 285 "csytri_3x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 286 "csytri_3x.f"
		i__1 = k + invd * work_dim1;
#line 286 "csytri_3x.f"
		z_div(&z__1, &akp1, &d__);
#line 286 "csytri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 287 "csytri_3x.f"
		i__1 = k + 1 + (invd + 1) * work_dim1;
#line 287 "csytri_3x.f"
		z_div(&z__1, &ak, &d__);
#line 287 "csytri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 288 "csytri_3x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 288 "csytri_3x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 288 "csytri_3x.f"
		z_div(&z__1, &z__2, &d__);
#line 288 "csytri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 289 "csytri_3x.f"
		i__1 = k + 1 + invd * work_dim1;
#line 289 "csytri_3x.f"
		i__2 = k + (invd + 1) * work_dim1;
#line 289 "csytri_3x.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 290 "csytri_3x.f"
		++k;
#line 291 "csytri_3x.f"
	    }
#line 292 "csytri_3x.f"
	    ++k;
#line 293 "csytri_3x.f"
	}

/*        inv(U**T) = (inv(U))**T */

/*        inv(U**T) * inv(D) * inv(U) */

#line 299 "csytri_3x.f"
	cut = *n;
#line 300 "csytri_3x.f"
	while(cut > 0) {
#line 301 "csytri_3x.f"
	    nnb = *nb;
#line 302 "csytri_3x.f"
	    if (cut <= nnb) {
#line 303 "csytri_3x.f"
		nnb = cut;
#line 304 "csytri_3x.f"
	    } else {
#line 305 "csytri_3x.f"
		icount = 0;
/*              count negative elements, */
#line 307 "csytri_3x.f"
		i__1 = cut;
#line 307 "csytri_3x.f"
		for (i__ = cut + 1 - nnb; i__ <= i__1; ++i__) {
#line 308 "csytri_3x.f"
		    if (ipiv[i__] < 0) {
#line 308 "csytri_3x.f"
			++icount;
#line 308 "csytri_3x.f"
		    }
#line 309 "csytri_3x.f"
		}
/*              need a even number for a clear cut */
#line 311 "csytri_3x.f"
		if (icount % 2 == 1) {
#line 311 "csytri_3x.f"
		    ++nnb;
#line 311 "csytri_3x.f"
		}
#line 312 "csytri_3x.f"
	    }
#line 314 "csytri_3x.f"
	    cut -= nnb;

/*           U01 Block */

#line 318 "csytri_3x.f"
	    i__1 = cut;
#line 318 "csytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 319 "csytri_3x.f"
		i__2 = nnb;
#line 319 "csytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 320 "csytri_3x.f"
		    i__3 = i__ + j * work_dim1;
#line 320 "csytri_3x.f"
		    i__4 = i__ + (cut + j) * a_dim1;
#line 320 "csytri_3x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 321 "csytri_3x.f"
		}
#line 322 "csytri_3x.f"
	    }

/*           U11 Block */

#line 326 "csytri_3x.f"
	    i__1 = nnb;
#line 326 "csytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 327 "csytri_3x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 327 "csytri_3x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 328 "csytri_3x.f"
		i__2 = i__ - 1;
#line 328 "csytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 329 "csytri_3x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 329 "csytri_3x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 330 "csytri_3x.f"
		}
#line 331 "csytri_3x.f"
		i__2 = nnb;
#line 331 "csytri_3x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 332 "csytri_3x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 332 "csytri_3x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 332 "csytri_3x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 333 "csytri_3x.f"
		}
#line 334 "csytri_3x.f"
	    }

/*           invD * U01 */

#line 338 "csytri_3x.f"
	    i__ = 1;
#line 339 "csytri_3x.f"
	    while(i__ <= cut) {
#line 340 "csytri_3x.f"
		if (ipiv[i__] > 0) {
#line 341 "csytri_3x.f"
		    i__1 = nnb;
#line 341 "csytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 342 "csytri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 342 "csytri_3x.f"
			i__3 = i__ + invd * work_dim1;
#line 342 "csytri_3x.f"
			i__4 = i__ + j * work_dim1;
#line 342 "csytri_3x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 342 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 343 "csytri_3x.f"
		    }
#line 344 "csytri_3x.f"
		} else {
#line 345 "csytri_3x.f"
		    i__1 = nnb;
#line 345 "csytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 346 "csytri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 346 "csytri_3x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 347 "csytri_3x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 347 "csytri_3x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 348 "csytri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 348 "csytri_3x.f"
			i__3 = i__ + invd * work_dim1;
#line 348 "csytri_3x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 348 "csytri_3x.f"
			i__4 = i__ + (invd + 1) * work_dim1;
#line 348 "csytri_3x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 348 "csytri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 348 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 350 "csytri_3x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 350 "csytri_3x.f"
			i__3 = i__ + 1 + invd * work_dim1;
#line 350 "csytri_3x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 350 "csytri_3x.f"
			i__4 = i__ + 1 + (invd + 1) * work_dim1;
#line 350 "csytri_3x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 350 "csytri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 350 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 352 "csytri_3x.f"
		    }
#line 353 "csytri_3x.f"
		    ++i__;
#line 354 "csytri_3x.f"
		}
#line 355 "csytri_3x.f"
		++i__;
#line 356 "csytri_3x.f"
	    }

/*           invD1 * U11 */

#line 360 "csytri_3x.f"
	    i__ = 1;
#line 361 "csytri_3x.f"
	    while(i__ <= nnb) {
#line 362 "csytri_3x.f"
		if (ipiv[cut + i__] > 0) {
#line 363 "csytri_3x.f"
		    i__1 = nnb;
#line 363 "csytri_3x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 364 "csytri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 364 "csytri_3x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 364 "csytri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 364 "csytri_3x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 364 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 365 "csytri_3x.f"
		    }
#line 366 "csytri_3x.f"
		} else {
#line 367 "csytri_3x.f"
		    i__1 = nnb;
#line 367 "csytri_3x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 368 "csytri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 368 "csytri_3x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 369 "csytri_3x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 369 "csytri_3x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 370 "csytri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 370 "csytri_3x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 370 "csytri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 370 "csytri_3x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 370 "csytri_3x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 370 "csytri_3x.f"
			i__6 = u11 + i__ + 1 + j * work_dim1;
#line 370 "csytri_3x.f"
			z__3.r = work[i__5].r * work[i__6].r - work[i__5].i * 
				work[i__6].i, z__3.i = work[i__5].r * work[
				i__6].i + work[i__5].i * work[i__6].r;
#line 370 "csytri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 370 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 372 "csytri_3x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 372 "csytri_3x.f"
			i__3 = cut + i__ + 1 + invd * work_dim1;
#line 372 "csytri_3x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 372 "csytri_3x.f"
			i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
#line 372 "csytri_3x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 372 "csytri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 372 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 374 "csytri_3x.f"
		    }
#line 375 "csytri_3x.f"
		    ++i__;
#line 376 "csytri_3x.f"
		}
#line 377 "csytri_3x.f"
		++i__;
#line 378 "csytri_3x.f"
	    }

/*           U11**T * invD1 * U11 -> U11 */

#line 382 "csytri_3x.f"
	    i__1 = *n + *nb + 1;
#line 382 "csytri_3x.f"
	    ctrmm_("L", "U", "T", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 
		    1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 386 "csytri_3x.f"
	    i__1 = nnb;
#line 386 "csytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 387 "csytri_3x.f"
		i__2 = nnb;
#line 387 "csytri_3x.f"
		for (j = i__; j <= i__2; ++j) {
#line 388 "csytri_3x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 388 "csytri_3x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 388 "csytri_3x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 389 "csytri_3x.f"
		}
#line 390 "csytri_3x.f"
	    }

/*           U01**T * invD * U01 -> A( CUT+I, CUT+J ) */

#line 394 "csytri_3x.f"
	    i__1 = *n + *nb + 1;
#line 394 "csytri_3x.f"
	    i__2 = *n + *nb + 1;
#line 394 "csytri_3x.f"
	    cgemm_("T", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 
		    1], lda, &work[work_offset], &i__1, &c_b2, &work[u11 + 1 
		    + work_dim1], &i__2, (ftnlen)1, (ftnlen)1);

/*           U11 =  U11**T * invD1 * U11 + U01**T * invD * U01 */

#line 401 "csytri_3x.f"
	    i__1 = nnb;
#line 401 "csytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 402 "csytri_3x.f"
		i__2 = nnb;
#line 402 "csytri_3x.f"
		for (j = i__; j <= i__2; ++j) {
#line 403 "csytri_3x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 403 "csytri_3x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 403 "csytri_3x.f"
		    i__5 = u11 + i__ + j * work_dim1;
#line 403 "csytri_3x.f"
		    z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i + 
			    work[i__5].i;
#line 403 "csytri_3x.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 404 "csytri_3x.f"
		}
#line 405 "csytri_3x.f"
	    }

/*           U01 =  U00**T * invD0 * U01 */

#line 409 "csytri_3x.f"
	    i__1 = *n + *nb + 1;
#line 409 "csytri_3x.f"
	    ctrmm_("L", uplo, "T", "U", &cut, &nnb, &c_b1, &a[a_offset], lda, 
		    &work[work_offset], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);

/*           Update U01 */

#line 415 "csytri_3x.f"
	    i__1 = cut;
#line 415 "csytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 416 "csytri_3x.f"
		i__2 = nnb;
#line 416 "csytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 417 "csytri_3x.f"
		    i__3 = i__ + (cut + j) * a_dim1;
#line 417 "csytri_3x.f"
		    i__4 = i__ + j * work_dim1;
#line 417 "csytri_3x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 418 "csytri_3x.f"
		}
#line 419 "csytri_3x.f"
	    }

/*           Next Block */

#line 423 "csytri_3x.f"
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

#line 436 "csytri_3x.f"
	i__1 = *n;
#line 436 "csytri_3x.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 437 "csytri_3x.f"
	    ip = (i__2 = ipiv[i__], abs(i__2));
#line 438 "csytri_3x.f"
	    if (ip != i__) {
#line 439 "csytri_3x.f"
		if (i__ < ip) {
#line 439 "csytri_3x.f"
		    csyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 439 "csytri_3x.f"
		}
#line 440 "csytri_3x.f"
		if (i__ > ip) {
#line 440 "csytri_3x.f"
		    csyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 440 "csytri_3x.f"
		}
#line 441 "csytri_3x.f"
	    }
#line 442 "csytri_3x.f"
	}

#line 444 "csytri_3x.f"
    } else {

/*        Begin Lower */

/*        inv A = P * inv(L**T) * inv(D) * inv(L) * P**T. */

#line 450 "csytri_3x.f"
	ctrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*        inv(D) and inv(D) * inv(L) */

#line 454 "csytri_3x.f"
	k = *n;
#line 455 "csytri_3x.f"
	while(k >= 1) {
#line 456 "csytri_3x.f"
	    if (ipiv[k] > 0) {
/*              1 x 1 diagonal NNB */
#line 458 "csytri_3x.f"
		i__1 = k + invd * work_dim1;
#line 458 "csytri_3x.f"
		z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 458 "csytri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 459 "csytri_3x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 459 "csytri_3x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 460 "csytri_3x.f"
	    } else {
/*              2 x 2 diagonal NNB */
#line 462 "csytri_3x.f"
		i__1 = k - 1 + work_dim1;
#line 462 "csytri_3x.f"
		t.r = work[i__1].r, t.i = work[i__1].i;
#line 463 "csytri_3x.f"
		z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &t);
#line 463 "csytri_3x.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 464 "csytri_3x.f"
		z_div(&z__1, &a[k + k * a_dim1], &t);
#line 464 "csytri_3x.f"
		akp1.r = z__1.r, akp1.i = z__1.i;
#line 465 "csytri_3x.f"
		z_div(&z__1, &work[k - 1 + work_dim1], &t);
#line 465 "csytri_3x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 466 "csytri_3x.f"
		z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * 
			akp1.i + ak.i * akp1.r;
#line 466 "csytri_3x.f"
		z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 466 "csytri_3x.f"
		z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + 
			t.i * z__2.r;
#line 466 "csytri_3x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 467 "csytri_3x.f"
		i__1 = k - 1 + invd * work_dim1;
#line 467 "csytri_3x.f"
		z_div(&z__1, &akp1, &d__);
#line 467 "csytri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 468 "csytri_3x.f"
		i__1 = k + invd * work_dim1;
#line 468 "csytri_3x.f"
		z_div(&z__1, &ak, &d__);
#line 468 "csytri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 469 "csytri_3x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 469 "csytri_3x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 469 "csytri_3x.f"
		z_div(&z__1, &z__2, &d__);
#line 469 "csytri_3x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 470 "csytri_3x.f"
		i__1 = k - 1 + (invd + 1) * work_dim1;
#line 470 "csytri_3x.f"
		i__2 = k + (invd + 1) * work_dim1;
#line 470 "csytri_3x.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 471 "csytri_3x.f"
		--k;
#line 472 "csytri_3x.f"
	    }
#line 473 "csytri_3x.f"
	    --k;
#line 474 "csytri_3x.f"
	}

/*        inv(L**T) = (inv(L))**T */

/*        inv(L**T) * inv(D) * inv(L) */

#line 480 "csytri_3x.f"
	cut = 0;
#line 481 "csytri_3x.f"
	while(cut < *n) {
#line 482 "csytri_3x.f"
	    nnb = *nb;
#line 483 "csytri_3x.f"
	    if (cut + nnb > *n) {
#line 484 "csytri_3x.f"
		nnb = *n - cut;
#line 485 "csytri_3x.f"
	    } else {
#line 486 "csytri_3x.f"
		icount = 0;
/*              count negative elements, */
#line 488 "csytri_3x.f"
		i__1 = cut + nnb;
#line 488 "csytri_3x.f"
		for (i__ = cut + 1; i__ <= i__1; ++i__) {
#line 489 "csytri_3x.f"
		    if (ipiv[i__] < 0) {
#line 489 "csytri_3x.f"
			++icount;
#line 489 "csytri_3x.f"
		    }
#line 490 "csytri_3x.f"
		}
/*              need a even number for a clear cut */
#line 492 "csytri_3x.f"
		if (icount % 2 == 1) {
#line 492 "csytri_3x.f"
		    ++nnb;
#line 492 "csytri_3x.f"
		}
#line 493 "csytri_3x.f"
	    }

/*           L21 Block */

#line 497 "csytri_3x.f"
	    i__1 = *n - cut - nnb;
#line 497 "csytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 498 "csytri_3x.f"
		i__2 = nnb;
#line 498 "csytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 499 "csytri_3x.f"
		    i__3 = i__ + j * work_dim1;
#line 499 "csytri_3x.f"
		    i__4 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 499 "csytri_3x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 500 "csytri_3x.f"
		}
#line 501 "csytri_3x.f"
	    }

/*           L11 Block */

#line 505 "csytri_3x.f"
	    i__1 = nnb;
#line 505 "csytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 506 "csytri_3x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 506 "csytri_3x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 507 "csytri_3x.f"
		i__2 = nnb;
#line 507 "csytri_3x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 508 "csytri_3x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 508 "csytri_3x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 509 "csytri_3x.f"
		}
#line 510 "csytri_3x.f"
		i__2 = i__ - 1;
#line 510 "csytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 511 "csytri_3x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 511 "csytri_3x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 511 "csytri_3x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 512 "csytri_3x.f"
		}
#line 513 "csytri_3x.f"
	    }

/*           invD*L21 */

#line 517 "csytri_3x.f"
	    i__ = *n - cut - nnb;
#line 518 "csytri_3x.f"
	    while(i__ >= 1) {
#line 519 "csytri_3x.f"
		if (ipiv[cut + nnb + i__] > 0) {
#line 520 "csytri_3x.f"
		    i__1 = nnb;
#line 520 "csytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 521 "csytri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 521 "csytri_3x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 521 "csytri_3x.f"
			i__4 = i__ + j * work_dim1;
#line 521 "csytri_3x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 521 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 522 "csytri_3x.f"
		    }
#line 523 "csytri_3x.f"
		} else {
#line 524 "csytri_3x.f"
		    i__1 = nnb;
#line 524 "csytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 525 "csytri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 525 "csytri_3x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 526 "csytri_3x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 526 "csytri_3x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 527 "csytri_3x.f"
			i__2 = i__ + j * work_dim1;
#line 527 "csytri_3x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 527 "csytri_3x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 527 "csytri_3x.f"
			i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
#line 527 "csytri_3x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 527 "csytri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 527 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 529 "csytri_3x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 529 "csytri_3x.f"
			i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
#line 529 "csytri_3x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 529 "csytri_3x.f"
			i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
#line 529 "csytri_3x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 529 "csytri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 529 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 531 "csytri_3x.f"
		    }
#line 532 "csytri_3x.f"
		    --i__;
#line 533 "csytri_3x.f"
		}
#line 534 "csytri_3x.f"
		--i__;
#line 535 "csytri_3x.f"
	    }

/*           invD1*L11 */

#line 539 "csytri_3x.f"
	    i__ = nnb;
#line 540 "csytri_3x.f"
	    while(i__ >= 1) {
#line 541 "csytri_3x.f"
		if (ipiv[cut + i__] > 0) {
#line 542 "csytri_3x.f"
		    i__1 = nnb;
#line 542 "csytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 543 "csytri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 543 "csytri_3x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 543 "csytri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 543 "csytri_3x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 543 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 544 "csytri_3x.f"
		    }
#line 546 "csytri_3x.f"
		} else {
#line 547 "csytri_3x.f"
		    i__1 = nnb;
#line 547 "csytri_3x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 548 "csytri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 548 "csytri_3x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 549 "csytri_3x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 549 "csytri_3x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 550 "csytri_3x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 550 "csytri_3x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 550 "csytri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 550 "csytri_3x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 550 "csytri_3x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 550 "csytri_3x.f"
			z__3.r = work[i__5].r * u11_ip1_j__.r - work[i__5].i *
				 u11_ip1_j__.i, z__3.i = work[i__5].r * 
				u11_ip1_j__.i + work[i__5].i * u11_ip1_j__.r;
#line 550 "csytri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 550 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 552 "csytri_3x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 552 "csytri_3x.f"
			i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
#line 552 "csytri_3x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 552 "csytri_3x.f"
			i__4 = cut + i__ - 1 + invd * work_dim1;
#line 552 "csytri_3x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 552 "csytri_3x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 552 "csytri_3x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 554 "csytri_3x.f"
		    }
#line 555 "csytri_3x.f"
		    --i__;
#line 556 "csytri_3x.f"
		}
#line 557 "csytri_3x.f"
		--i__;
#line 558 "csytri_3x.f"
	    }

/*           L11**T * invD1 * L11 -> L11 */

#line 562 "csytri_3x.f"
	    i__1 = *n + *nb + 1;
#line 562 "csytri_3x.f"
	    ctrmm_("L", uplo, "T", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 567 "csytri_3x.f"
	    i__1 = nnb;
#line 567 "csytri_3x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 568 "csytri_3x.f"
		i__2 = i__;
#line 568 "csytri_3x.f"
		for (j = 1; j <= i__2; ++j) {
#line 569 "csytri_3x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 569 "csytri_3x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 569 "csytri_3x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 570 "csytri_3x.f"
		}
#line 571 "csytri_3x.f"
	    }

#line 573 "csytri_3x.f"
	    if (cut + nnb < *n) {

/*              L21**T * invD2*L21 -> A( CUT+I, CUT+J ) */

#line 577 "csytri_3x.f"
		i__1 = *n - nnb - cut;
#line 577 "csytri_3x.f"
		i__2 = *n + *nb + 1;
#line 577 "csytri_3x.f"
		i__3 = *n + *nb + 1;
#line 577 "csytri_3x.f"
		cgemm_("T", "N", &nnb, &nnb, &i__1, &c_b1, &a[cut + nnb + 1 + 
			(cut + 1) * a_dim1], lda, &work[work_offset], &i__2, &
			c_b2, &work[u11 + 1 + work_dim1], &i__3, (ftnlen)1, (
			ftnlen)1);

/*              L11 =  L11**T * invD1 * L11 + U01**T * invD * U01 */

#line 584 "csytri_3x.f"
		i__1 = nnb;
#line 584 "csytri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 585 "csytri_3x.f"
		    i__2 = i__;
#line 585 "csytri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 586 "csytri_3x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 586 "csytri_3x.f"
			i__4 = cut + i__ + (cut + j) * a_dim1;
#line 586 "csytri_3x.f"
			i__5 = u11 + i__ + j * work_dim1;
#line 586 "csytri_3x.f"
			z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i 
				+ work[i__5].i;
#line 586 "csytri_3x.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 587 "csytri_3x.f"
		    }
#line 588 "csytri_3x.f"
		}

/*              L01 =  L22**T * invD2 * L21 */

#line 592 "csytri_3x.f"
		i__1 = *n - nnb - cut;
#line 592 "csytri_3x.f"
		i__2 = *n + *nb + 1;
#line 592 "csytri_3x.f"
		ctrmm_("L", uplo, "T", "U", &i__1, &nnb, &c_b1, &a[cut + nnb 
			+ 1 + (cut + nnb + 1) * a_dim1], lda, &work[
			work_offset], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);

/*              Update L21 */

#line 598 "csytri_3x.f"
		i__1 = *n - cut - nnb;
#line 598 "csytri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 599 "csytri_3x.f"
		    i__2 = nnb;
#line 599 "csytri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 600 "csytri_3x.f"
			i__3 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 600 "csytri_3x.f"
			i__4 = i__ + j * work_dim1;
#line 600 "csytri_3x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 601 "csytri_3x.f"
		    }
#line 602 "csytri_3x.f"
		}

#line 604 "csytri_3x.f"
	    } else {

/*              L11 =  L11**T * invD1 * L11 */

#line 608 "csytri_3x.f"
		i__1 = nnb;
#line 608 "csytri_3x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 609 "csytri_3x.f"
		    i__2 = i__;
#line 609 "csytri_3x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 610 "csytri_3x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 610 "csytri_3x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 610 "csytri_3x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 611 "csytri_3x.f"
		    }
#line 612 "csytri_3x.f"
		}
#line 613 "csytri_3x.f"
	    }

/*           Next Block */

#line 617 "csytri_3x.f"
	    cut += nnb;

#line 619 "csytri_3x.f"
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

#line 632 "csytri_3x.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 633 "csytri_3x.f"
	    ip = (i__1 = ipiv[i__], abs(i__1));
#line 634 "csytri_3x.f"
	    if (ip != i__) {
#line 635 "csytri_3x.f"
		if (i__ < ip) {
#line 635 "csytri_3x.f"
		    csyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 635 "csytri_3x.f"
		}
#line 636 "csytri_3x.f"
		if (i__ > ip) {
#line 636 "csytri_3x.f"
		    csyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 636 "csytri_3x.f"
		}
#line 637 "csytri_3x.f"
	    }
#line 638 "csytri_3x.f"
	}

#line 640 "csytri_3x.f"
    }

#line 642 "csytri_3x.f"
    return 0;

/*     End of CSYTRI_3X */

} /* csytri_3x__ */

