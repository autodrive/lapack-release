#line 1 "zsytri2x.f"
/* zsytri2x.f -- translated by f2c (version 20100827).
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

#line 1 "zsytri2x.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};

/* > \brief \b ZSYTRI2X */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRI2X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri2
x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri2
x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri2
x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( N+NB+1,* ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYTRI2X computes the inverse of a complex symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > ZSYTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
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
/* >          On entry, the NNB diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by ZSYTRF. */
/* > */
/* >          On exit, if INFO = 0, the (symmetric) inverse of the original */
/* >          matrix.  If UPLO = 'U', the upper triangular part of the */
/* >          inverse is formed and the part of A below the diagonal is not */
/* >          referenced; if UPLO = 'L' the lower triangular part of the */
/* >          inverse is formed and the part of A above the diagonal is */
/* >          not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the NNB structure of D */
/* >          as determined by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N+NNB+1,NNB+3) */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          Block size */
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

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsytri2x_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *nb, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex d__;
    static integer i__, j, k;
    extern /* Subroutine */ int zsyswapr_(char *, integer *, doublecomplex *, 
	    integer *, integer *, integer *, ftnlen);
    static doublecomplex t, ak;
    static integer u11, ip, nnb, cut;
    static doublecomplex akp1;
    static integer invd;
    static doublecomplex akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer count;
    static logical upper;
    extern /* Subroutine */ int ztrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublecomplex u01_i_j__, u11_i_j__;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), ztrtri_(
	    char *, char *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen, ftnlen);
    static doublecomplex u01_ip1_j__, u11_ip1_j__;
    extern /* Subroutine */ int zsyconv_(char *, char *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     ftnlen, ftnlen);


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

#line 169 "zsytri2x.f"
    /* Parameter adjustments */
#line 169 "zsytri2x.f"
    a_dim1 = *lda;
#line 169 "zsytri2x.f"
    a_offset = 1 + a_dim1;
#line 169 "zsytri2x.f"
    a -= a_offset;
#line 169 "zsytri2x.f"
    --ipiv;
#line 169 "zsytri2x.f"
    work_dim1 = *n + *nb + 1;
#line 169 "zsytri2x.f"
    work_offset = 1 + work_dim1;
#line 169 "zsytri2x.f"
    work -= work_offset;
#line 169 "zsytri2x.f"

#line 169 "zsytri2x.f"
    /* Function Body */
#line 169 "zsytri2x.f"
    *info = 0;
#line 170 "zsytri2x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 171 "zsytri2x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 172 "zsytri2x.f"
	*info = -1;
#line 173 "zsytri2x.f"
    } else if (*n < 0) {
#line 174 "zsytri2x.f"
	*info = -2;
#line 175 "zsytri2x.f"
    } else if (*lda < max(1,*n)) {
#line 176 "zsytri2x.f"
	*info = -4;
#line 177 "zsytri2x.f"
    }

/*     Quick return if possible */


#line 182 "zsytri2x.f"
    if (*info != 0) {
#line 183 "zsytri2x.f"
	i__1 = -(*info);
#line 183 "zsytri2x.f"
	xerbla_("ZSYTRI2X", &i__1, (ftnlen)8);
#line 184 "zsytri2x.f"
	return 0;
#line 185 "zsytri2x.f"
    }
#line 186 "zsytri2x.f"
    if (*n == 0) {
#line 186 "zsytri2x.f"
	return 0;
#line 186 "zsytri2x.f"
    }

/*     Convert A */
/*     Workspace got Non-diag elements of D */

#line 192 "zsytri2x.f"
    zsyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     Check that the diagonal matrix D is nonsingular. */

#line 196 "zsytri2x.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 200 "zsytri2x.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 201 "zsytri2x.f"
	    i__1 = *info + *info * a_dim1;
#line 201 "zsytri2x.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 201 "zsytri2x.f"
		return 0;
#line 201 "zsytri2x.f"
	    }
#line 203 "zsytri2x.f"
	}
#line 204 "zsytri2x.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 208 "zsytri2x.f"
	i__1 = *n;
#line 208 "zsytri2x.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 209 "zsytri2x.f"
	    i__2 = *info + *info * a_dim1;
#line 209 "zsytri2x.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 209 "zsytri2x.f"
		return 0;
#line 209 "zsytri2x.f"
	    }
#line 211 "zsytri2x.f"
	}
#line 212 "zsytri2x.f"
    }
#line 213 "zsytri2x.f"
    *info = 0;

/*  Splitting Workspace */
/*     U01 is a block (N,NB+1) */
/*     The first element of U01 is in WORK(1,1) */
/*     U11 is a block (NB+1,NB+1) */
/*     The first element of U11 is in WORK(N+1,1) */
#line 220 "zsytri2x.f"
    u11 = *n;
/*     INVD is a block (N,2) */
/*     The first element of INVD is in WORK(1,INVD) */
#line 223 "zsytri2x.f"
    invd = *nb + 2;
#line 225 "zsytri2x.f"
    if (upper) {

/*        invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */

#line 229 "zsytri2x.f"
	ztrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 233 "zsytri2x.f"
	k = 1;
#line 234 "zsytri2x.f"
	while(k <= *n) {
#line 235 "zsytri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 237 "zsytri2x.f"
		i__1 = k + invd * work_dim1;
#line 237 "zsytri2x.f"
		z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 237 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 238 "zsytri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 238 "zsytri2x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 239 "zsytri2x.f"
		++k;
#line 240 "zsytri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 242 "zsytri2x.f"
		i__1 = k + 1 + work_dim1;
#line 242 "zsytri2x.f"
		t.r = work[i__1].r, t.i = work[i__1].i;
#line 243 "zsytri2x.f"
		z_div(&z__1, &a[k + k * a_dim1], &t);
#line 243 "zsytri2x.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 244 "zsytri2x.f"
		z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &t);
#line 244 "zsytri2x.f"
		akp1.r = z__1.r, akp1.i = z__1.i;
#line 245 "zsytri2x.f"
		z_div(&z__1, &work[k + 1 + work_dim1], &t);
#line 245 "zsytri2x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 246 "zsytri2x.f"
		z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * 
			akp1.i + ak.i * akp1.r;
#line 246 "zsytri2x.f"
		z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 246 "zsytri2x.f"
		z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + 
			t.i * z__2.r;
#line 246 "zsytri2x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 247 "zsytri2x.f"
		i__1 = k + invd * work_dim1;
#line 247 "zsytri2x.f"
		z_div(&z__1, &akp1, &d__);
#line 247 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 248 "zsytri2x.f"
		i__1 = k + 1 + (invd + 1) * work_dim1;
#line 248 "zsytri2x.f"
		z_div(&z__1, &ak, &d__);
#line 248 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 249 "zsytri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 249 "zsytri2x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 249 "zsytri2x.f"
		z_div(&z__1, &z__2, &d__);
#line 249 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 250 "zsytri2x.f"
		i__1 = k + 1 + invd * work_dim1;
#line 250 "zsytri2x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 250 "zsytri2x.f"
		z_div(&z__1, &z__2, &d__);
#line 250 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 251 "zsytri2x.f"
		k += 2;
#line 252 "zsytri2x.f"
	    }
#line 253 "zsytri2x.f"
	}

/*       inv(U**T) = (inv(U))**T */

/*       inv(U**T)*inv(D)*inv(U) */

#line 259 "zsytri2x.f"
	cut = *n;
#line 260 "zsytri2x.f"
	while(cut > 0) {
#line 261 "zsytri2x.f"
	    nnb = *nb;
#line 262 "zsytri2x.f"
	    if (cut <= nnb) {
#line 263 "zsytri2x.f"
		nnb = cut;
#line 264 "zsytri2x.f"
	    } else {
#line 265 "zsytri2x.f"
		count = 0;
/*             count negative elements, */
#line 267 "zsytri2x.f"
		i__1 = cut;
#line 267 "zsytri2x.f"
		for (i__ = cut + 1 - nnb; i__ <= i__1; ++i__) {
#line 268 "zsytri2x.f"
		    if (ipiv[i__] < 0) {
#line 268 "zsytri2x.f"
			++count;
#line 268 "zsytri2x.f"
		    }
#line 269 "zsytri2x.f"
		}
/*             need a even number for a clear cut */
#line 271 "zsytri2x.f"
		if (count % 2 == 1) {
#line 271 "zsytri2x.f"
		    ++nnb;
#line 271 "zsytri2x.f"
		}
#line 272 "zsytri2x.f"
	    }
#line 274 "zsytri2x.f"
	    cut -= nnb;

/*          U01 Block */

#line 278 "zsytri2x.f"
	    i__1 = cut;
#line 278 "zsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 279 "zsytri2x.f"
		i__2 = nnb;
#line 279 "zsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 280 "zsytri2x.f"
		    i__3 = i__ + j * work_dim1;
#line 280 "zsytri2x.f"
		    i__4 = i__ + (cut + j) * a_dim1;
#line 280 "zsytri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 281 "zsytri2x.f"
		}
#line 282 "zsytri2x.f"
	    }

/*          U11 Block */

#line 286 "zsytri2x.f"
	    i__1 = nnb;
#line 286 "zsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 287 "zsytri2x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 287 "zsytri2x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 288 "zsytri2x.f"
		i__2 = i__ - 1;
#line 288 "zsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 289 "zsytri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 289 "zsytri2x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 290 "zsytri2x.f"
		}
#line 291 "zsytri2x.f"
		i__2 = nnb;
#line 291 "zsytri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 292 "zsytri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 292 "zsytri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 292 "zsytri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 293 "zsytri2x.f"
		}
#line 294 "zsytri2x.f"
	    }

/*          invD*U01 */

#line 298 "zsytri2x.f"
	    i__ = 1;
#line 299 "zsytri2x.f"
	    while(i__ <= cut) {
#line 300 "zsytri2x.f"
		if (ipiv[i__] > 0) {
#line 301 "zsytri2x.f"
		    i__1 = nnb;
#line 301 "zsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 302 "zsytri2x.f"
			i__2 = i__ + j * work_dim1;
#line 302 "zsytri2x.f"
			i__3 = i__ + invd * work_dim1;
#line 302 "zsytri2x.f"
			i__4 = i__ + j * work_dim1;
#line 302 "zsytri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 302 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 303 "zsytri2x.f"
		    }
#line 304 "zsytri2x.f"
		    ++i__;
#line 305 "zsytri2x.f"
		} else {
#line 306 "zsytri2x.f"
		    i__1 = nnb;
#line 306 "zsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 307 "zsytri2x.f"
			i__2 = i__ + j * work_dim1;
#line 307 "zsytri2x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 308 "zsytri2x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 308 "zsytri2x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 309 "zsytri2x.f"
			i__2 = i__ + j * work_dim1;
#line 309 "zsytri2x.f"
			i__3 = i__ + invd * work_dim1;
#line 309 "zsytri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 309 "zsytri2x.f"
			i__4 = i__ + (invd + 1) * work_dim1;
#line 309 "zsytri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 309 "zsytri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 309 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 311 "zsytri2x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 311 "zsytri2x.f"
			i__3 = i__ + 1 + invd * work_dim1;
#line 311 "zsytri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 311 "zsytri2x.f"
			i__4 = i__ + 1 + (invd + 1) * work_dim1;
#line 311 "zsytri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 311 "zsytri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 311 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 313 "zsytri2x.f"
		    }
#line 314 "zsytri2x.f"
		    i__ += 2;
#line 315 "zsytri2x.f"
		}
#line 316 "zsytri2x.f"
	    }

/*        invD1*U11 */

#line 320 "zsytri2x.f"
	    i__ = 1;
#line 321 "zsytri2x.f"
	    while(i__ <= nnb) {
#line 322 "zsytri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 323 "zsytri2x.f"
		    i__1 = nnb;
#line 323 "zsytri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 324 "zsytri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 324 "zsytri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 324 "zsytri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 324 "zsytri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 324 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 325 "zsytri2x.f"
		    }
#line 326 "zsytri2x.f"
		    ++i__;
#line 327 "zsytri2x.f"
		} else {
#line 328 "zsytri2x.f"
		    i__1 = nnb;
#line 328 "zsytri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 329 "zsytri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 329 "zsytri2x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 330 "zsytri2x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 330 "zsytri2x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 331 "zsytri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 331 "zsytri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 331 "zsytri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 331 "zsytri2x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 331 "zsytri2x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 331 "zsytri2x.f"
			i__6 = u11 + i__ + 1 + j * work_dim1;
#line 331 "zsytri2x.f"
			z__3.r = work[i__5].r * work[i__6].r - work[i__5].i * 
				work[i__6].i, z__3.i = work[i__5].r * work[
				i__6].i + work[i__5].i * work[i__6].r;
#line 331 "zsytri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 331 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 333 "zsytri2x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 333 "zsytri2x.f"
			i__3 = cut + i__ + 1 + invd * work_dim1;
#line 333 "zsytri2x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 333 "zsytri2x.f"
			i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
#line 333 "zsytri2x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 333 "zsytri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 333 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 335 "zsytri2x.f"
		    }
#line 336 "zsytri2x.f"
		    i__ += 2;
#line 337 "zsytri2x.f"
		}
#line 338 "zsytri2x.f"
	    }

/*       U11**T*invD1*U11->U11 */

#line 342 "zsytri2x.f"
	    i__1 = *n + *nb + 1;
#line 342 "zsytri2x.f"
	    ztrmm_("L", "U", "T", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 
		    1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 345 "zsytri2x.f"
	    i__1 = nnb;
#line 345 "zsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 346 "zsytri2x.f"
		i__2 = nnb;
#line 346 "zsytri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 347 "zsytri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 347 "zsytri2x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 347 "zsytri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 348 "zsytri2x.f"
		}
#line 349 "zsytri2x.f"
	    }

/*          U01**T*invD*U01->A(CUT+I,CUT+J) */

#line 353 "zsytri2x.f"
	    i__1 = *n + *nb + 1;
#line 353 "zsytri2x.f"
	    i__2 = *n + *nb + 1;
#line 353 "zsytri2x.f"
	    zgemm_("T", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 
		    1], lda, &work[work_offset], &i__1, &c_b2, &work[u11 + 1 
		    + work_dim1], &i__2, (ftnlen)1, (ftnlen)1);

/*        U11 =  U11**T*invD1*U11 + U01**T*invD*U01 */

#line 358 "zsytri2x.f"
	    i__1 = nnb;
#line 358 "zsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 359 "zsytri2x.f"
		i__2 = nnb;
#line 359 "zsytri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 360 "zsytri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 360 "zsytri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 360 "zsytri2x.f"
		    i__5 = u11 + i__ + j * work_dim1;
#line 360 "zsytri2x.f"
		    z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i + 
			    work[i__5].i;
#line 360 "zsytri2x.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 361 "zsytri2x.f"
		}
#line 362 "zsytri2x.f"
	    }

/*        U01 =  U00**T*invD0*U01 */

#line 366 "zsytri2x.f"
	    i__1 = *n + *nb + 1;
#line 366 "zsytri2x.f"
	    ztrmm_("L", uplo, "T", "U", &cut, &nnb, &c_b1, &a[a_offset], lda, 
		    &work[work_offset], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);

/*        Update U01 */

#line 372 "zsytri2x.f"
	    i__1 = cut;
#line 372 "zsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 373 "zsytri2x.f"
		i__2 = nnb;
#line 373 "zsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 374 "zsytri2x.f"
		    i__3 = i__ + (cut + j) * a_dim1;
#line 374 "zsytri2x.f"
		    i__4 = i__ + j * work_dim1;
#line 374 "zsytri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 375 "zsytri2x.f"
		}
#line 376 "zsytri2x.f"
	    }

/*      Next Block */

#line 380 "zsytri2x.f"
	}

/*        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */

#line 384 "zsytri2x.f"
	i__ = 1;
#line 385 "zsytri2x.f"
	while(i__ <= *n) {
#line 386 "zsytri2x.f"
	    if (ipiv[i__] > 0) {
#line 387 "zsytri2x.f"
		ip = ipiv[i__];
#line 388 "zsytri2x.f"
		if (i__ < ip) {
#line 388 "zsytri2x.f"
		    zsyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 388 "zsytri2x.f"
		}
#line 389 "zsytri2x.f"
		if (i__ > ip) {
#line 389 "zsytri2x.f"
		    zsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 389 "zsytri2x.f"
		}
#line 390 "zsytri2x.f"
	    } else {
#line 391 "zsytri2x.f"
		ip = -ipiv[i__];
#line 392 "zsytri2x.f"
		++i__;
#line 393 "zsytri2x.f"
		if (i__ - 1 < ip) {
#line 393 "zsytri2x.f"
		    i__1 = i__ - 1;
#line 393 "zsytri2x.f"
		    zsyswapr_(uplo, n, &a[a_offset], lda, &i__1, &ip, (ftnlen)
			    1);
#line 393 "zsytri2x.f"
		}
#line 395 "zsytri2x.f"
		if (i__ - 1 > ip) {
#line 395 "zsytri2x.f"
		    i__1 = i__ - 1;
#line 395 "zsytri2x.f"
		    zsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__1, (ftnlen)
			    1);
#line 395 "zsytri2x.f"
		}
#line 397 "zsytri2x.f"
	    }
#line 398 "zsytri2x.f"
	    ++i__;
#line 399 "zsytri2x.f"
	}
#line 400 "zsytri2x.f"
    } else {

/*        LOWER... */

/*        invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */

#line 406 "zsytri2x.f"
	ztrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 410 "zsytri2x.f"
	k = *n;
#line 411 "zsytri2x.f"
	while(k >= 1) {
#line 412 "zsytri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 414 "zsytri2x.f"
		i__1 = k + invd * work_dim1;
#line 414 "zsytri2x.f"
		z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 414 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 415 "zsytri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 415 "zsytri2x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 416 "zsytri2x.f"
		--k;
#line 417 "zsytri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 419 "zsytri2x.f"
		i__1 = k - 1 + work_dim1;
#line 419 "zsytri2x.f"
		t.r = work[i__1].r, t.i = work[i__1].i;
#line 420 "zsytri2x.f"
		z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &t);
#line 420 "zsytri2x.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 421 "zsytri2x.f"
		z_div(&z__1, &a[k + k * a_dim1], &t);
#line 421 "zsytri2x.f"
		akp1.r = z__1.r, akp1.i = z__1.i;
#line 422 "zsytri2x.f"
		z_div(&z__1, &work[k - 1 + work_dim1], &t);
#line 422 "zsytri2x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 423 "zsytri2x.f"
		z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * 
			akp1.i + ak.i * akp1.r;
#line 423 "zsytri2x.f"
		z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 423 "zsytri2x.f"
		z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + 
			t.i * z__2.r;
#line 423 "zsytri2x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 424 "zsytri2x.f"
		i__1 = k - 1 + invd * work_dim1;
#line 424 "zsytri2x.f"
		z_div(&z__1, &akp1, &d__);
#line 424 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 425 "zsytri2x.f"
		i__1 = k + invd * work_dim1;
#line 425 "zsytri2x.f"
		z_div(&z__1, &ak, &d__);
#line 425 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 426 "zsytri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 426 "zsytri2x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 426 "zsytri2x.f"
		z_div(&z__1, &z__2, &d__);
#line 426 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 427 "zsytri2x.f"
		i__1 = k - 1 + (invd + 1) * work_dim1;
#line 427 "zsytri2x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 427 "zsytri2x.f"
		z_div(&z__1, &z__2, &d__);
#line 427 "zsytri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 428 "zsytri2x.f"
		k += -2;
#line 429 "zsytri2x.f"
	    }
#line 430 "zsytri2x.f"
	}

/*       inv(U**T) = (inv(U))**T */

/*       inv(U**T)*inv(D)*inv(U) */

#line 436 "zsytri2x.f"
	cut = 0;
#line 437 "zsytri2x.f"
	while(cut < *n) {
#line 438 "zsytri2x.f"
	    nnb = *nb;
#line 439 "zsytri2x.f"
	    if (cut + nnb >= *n) {
#line 440 "zsytri2x.f"
		nnb = *n - cut;
#line 441 "zsytri2x.f"
	    } else {
#line 442 "zsytri2x.f"
		count = 0;
/*             count negative elements, */
#line 444 "zsytri2x.f"
		i__1 = cut + nnb;
#line 444 "zsytri2x.f"
		for (i__ = cut + 1; i__ <= i__1; ++i__) {
#line 445 "zsytri2x.f"
		    if (ipiv[i__] < 0) {
#line 445 "zsytri2x.f"
			++count;
#line 445 "zsytri2x.f"
		    }
#line 446 "zsytri2x.f"
		}
/*             need a even number for a clear cut */
#line 448 "zsytri2x.f"
		if (count % 2 == 1) {
#line 448 "zsytri2x.f"
		    ++nnb;
#line 448 "zsytri2x.f"
		}
#line 449 "zsytri2x.f"
	    }
/*      L21 Block */
#line 451 "zsytri2x.f"
	    i__1 = *n - cut - nnb;
#line 451 "zsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 452 "zsytri2x.f"
		i__2 = nnb;
#line 452 "zsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 453 "zsytri2x.f"
		    i__3 = i__ + j * work_dim1;
#line 453 "zsytri2x.f"
		    i__4 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 453 "zsytri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 454 "zsytri2x.f"
		}
#line 455 "zsytri2x.f"
	    }
/*     L11 Block */
#line 457 "zsytri2x.f"
	    i__1 = nnb;
#line 457 "zsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 458 "zsytri2x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 458 "zsytri2x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 459 "zsytri2x.f"
		i__2 = nnb;
#line 459 "zsytri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 460 "zsytri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 460 "zsytri2x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 461 "zsytri2x.f"
		}
#line 462 "zsytri2x.f"
		i__2 = i__ - 1;
#line 462 "zsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 463 "zsytri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 463 "zsytri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 463 "zsytri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 464 "zsytri2x.f"
		}
#line 465 "zsytri2x.f"
	    }

/*          invD*L21 */

#line 469 "zsytri2x.f"
	    i__ = *n - cut - nnb;
#line 470 "zsytri2x.f"
	    while(i__ >= 1) {
#line 471 "zsytri2x.f"
		if (ipiv[cut + nnb + i__] > 0) {
#line 472 "zsytri2x.f"
		    i__1 = nnb;
#line 472 "zsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 473 "zsytri2x.f"
			i__2 = i__ + j * work_dim1;
#line 473 "zsytri2x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 473 "zsytri2x.f"
			i__4 = i__ + j * work_dim1;
#line 473 "zsytri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 473 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 474 "zsytri2x.f"
		    }
#line 475 "zsytri2x.f"
		    --i__;
#line 476 "zsytri2x.f"
		} else {
#line 477 "zsytri2x.f"
		    i__1 = nnb;
#line 477 "zsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 478 "zsytri2x.f"
			i__2 = i__ + j * work_dim1;
#line 478 "zsytri2x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 479 "zsytri2x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 479 "zsytri2x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 480 "zsytri2x.f"
			i__2 = i__ + j * work_dim1;
#line 480 "zsytri2x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 480 "zsytri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 480 "zsytri2x.f"
			i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
#line 480 "zsytri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 480 "zsytri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 480 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 482 "zsytri2x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 482 "zsytri2x.f"
			i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
#line 482 "zsytri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 482 "zsytri2x.f"
			i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
#line 482 "zsytri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 482 "zsytri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 482 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 484 "zsytri2x.f"
		    }
#line 485 "zsytri2x.f"
		    i__ += -2;
#line 486 "zsytri2x.f"
		}
#line 487 "zsytri2x.f"
	    }

/*        invD1*L11 */

#line 491 "zsytri2x.f"
	    i__ = nnb;
#line 492 "zsytri2x.f"
	    while(i__ >= 1) {
#line 493 "zsytri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 494 "zsytri2x.f"
		    i__1 = nnb;
#line 494 "zsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 495 "zsytri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 495 "zsytri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 495 "zsytri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 495 "zsytri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 495 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 496 "zsytri2x.f"
		    }
#line 497 "zsytri2x.f"
		    --i__;
#line 498 "zsytri2x.f"
		} else {
#line 499 "zsytri2x.f"
		    i__1 = nnb;
#line 499 "zsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 500 "zsytri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 500 "zsytri2x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 501 "zsytri2x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 501 "zsytri2x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 502 "zsytri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 502 "zsytri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 502 "zsytri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 502 "zsytri2x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 502 "zsytri2x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 502 "zsytri2x.f"
			z__3.r = work[i__5].r * u11_ip1_j__.r - work[i__5].i *
				 u11_ip1_j__.i, z__3.i = work[i__5].r * 
				u11_ip1_j__.i + work[i__5].i * u11_ip1_j__.r;
#line 502 "zsytri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 502 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 504 "zsytri2x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 504 "zsytri2x.f"
			i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
#line 504 "zsytri2x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 504 "zsytri2x.f"
			i__4 = cut + i__ - 1 + invd * work_dim1;
#line 504 "zsytri2x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 504 "zsytri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 504 "zsytri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 506 "zsytri2x.f"
		    }
#line 507 "zsytri2x.f"
		    i__ += -2;
#line 508 "zsytri2x.f"
		}
#line 509 "zsytri2x.f"
	    }

/*       L11**T*invD1*L11->L11 */

#line 513 "zsytri2x.f"
	    i__1 = *n + *nb + 1;
#line 513 "zsytri2x.f"
	    ztrmm_("L", uplo, "T", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 516 "zsytri2x.f"
	    i__1 = nnb;
#line 516 "zsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 517 "zsytri2x.f"
		i__2 = i__;
#line 517 "zsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 518 "zsytri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 518 "zsytri2x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 518 "zsytri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 519 "zsytri2x.f"
		}
#line 520 "zsytri2x.f"
	    }

#line 523 "zsytri2x.f"
	    if (cut + nnb < *n) {

/*          L21**T*invD2*L21->A(CUT+I,CUT+J) */

#line 527 "zsytri2x.f"
		i__1 = *n - nnb - cut;
#line 527 "zsytri2x.f"
		i__2 = *n + *nb + 1;
#line 527 "zsytri2x.f"
		i__3 = *n + *nb + 1;
#line 527 "zsytri2x.f"
		zgemm_("T", "N", &nnb, &nnb, &i__1, &c_b1, &a[cut + nnb + 1 + 
			(cut + 1) * a_dim1], lda, &work[work_offset], &i__2, &
			c_b2, &work[u11 + 1 + work_dim1], &i__3, (ftnlen)1, (
			ftnlen)1);

/*        L11 =  L11**T*invD1*L11 + U01**T*invD*U01 */

#line 533 "zsytri2x.f"
		i__1 = nnb;
#line 533 "zsytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 534 "zsytri2x.f"
		    i__2 = i__;
#line 534 "zsytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 535 "zsytri2x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 535 "zsytri2x.f"
			i__4 = cut + i__ + (cut + j) * a_dim1;
#line 535 "zsytri2x.f"
			i__5 = u11 + i__ + j * work_dim1;
#line 535 "zsytri2x.f"
			z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i 
				+ work[i__5].i;
#line 535 "zsytri2x.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 536 "zsytri2x.f"
		    }
#line 537 "zsytri2x.f"
		}

/*        U01 =  L22**T*invD2*L21 */

#line 541 "zsytri2x.f"
		i__1 = *n - nnb - cut;
#line 541 "zsytri2x.f"
		i__2 = *n + *nb + 1;
#line 541 "zsytri2x.f"
		ztrmm_("L", uplo, "T", "U", &i__1, &nnb, &c_b1, &a[cut + nnb 
			+ 1 + (cut + nnb + 1) * a_dim1], lda, &work[
			work_offset], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
/*      Update L21 */
#line 545 "zsytri2x.f"
		i__1 = *n - cut - nnb;
#line 545 "zsytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 546 "zsytri2x.f"
		    i__2 = nnb;
#line 546 "zsytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 547 "zsytri2x.f"
			i__3 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 547 "zsytri2x.f"
			i__4 = i__ + j * work_dim1;
#line 547 "zsytri2x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 548 "zsytri2x.f"
		    }
#line 549 "zsytri2x.f"
		}
#line 550 "zsytri2x.f"
	    } else {

/*        L11 =  L11**T*invD1*L11 */

#line 554 "zsytri2x.f"
		i__1 = nnb;
#line 554 "zsytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 555 "zsytri2x.f"
		    i__2 = i__;
#line 555 "zsytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 556 "zsytri2x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 556 "zsytri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 556 "zsytri2x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 557 "zsytri2x.f"
		    }
#line 558 "zsytri2x.f"
		}
#line 559 "zsytri2x.f"
	    }

/*      Next Block */

#line 563 "zsytri2x.f"
	    cut += nnb;
#line 564 "zsytri2x.f"
	}

/*        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */

#line 568 "zsytri2x.f"
	i__ = *n;
#line 569 "zsytri2x.f"
	while(i__ >= 1) {
#line 570 "zsytri2x.f"
	    if (ipiv[i__] > 0) {
#line 571 "zsytri2x.f"
		ip = ipiv[i__];
#line 572 "zsytri2x.f"
		if (i__ < ip) {
#line 572 "zsytri2x.f"
		    zsyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 572 "zsytri2x.f"
		}
#line 573 "zsytri2x.f"
		if (i__ > ip) {
#line 573 "zsytri2x.f"
		    zsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 573 "zsytri2x.f"
		}
#line 574 "zsytri2x.f"
	    } else {
#line 575 "zsytri2x.f"
		ip = -ipiv[i__];
#line 576 "zsytri2x.f"
		if (i__ < ip) {
#line 576 "zsytri2x.f"
		    zsyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 576 "zsytri2x.f"
		}
#line 577 "zsytri2x.f"
		if (i__ > ip) {
#line 577 "zsytri2x.f"
		    zsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 577 "zsytri2x.f"
		}
#line 578 "zsytri2x.f"
		--i__;
#line 579 "zsytri2x.f"
	    }
#line 580 "zsytri2x.f"
	    --i__;
#line 581 "zsytri2x.f"
	}
#line 582 "zsytri2x.f"
    }

#line 584 "zsytri2x.f"
    return 0;

/*     End of ZSYTRI2X */

} /* zsytri2x_ */

