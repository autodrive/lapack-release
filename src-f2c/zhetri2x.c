#line 1 "zhetri2x.f"
/* zhetri2x.f -- translated by f2c (version 20100827).
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

#line 1 "zhetri2x.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};

/* > \brief \b ZHETRI2X */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRI2X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri2
x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri2
x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri2
x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16            A( LDA, * ), WORK( N+NB+1,* ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRI2X computes the inverse of a COMPLEX*16 Hermitian indefinite matrix */
/* > A using the factorization A = U*D*U**H or A = L*D*L**H computed by */
/* > ZHETRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**H; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**H. */
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
/* >          used to obtain the factor U or L as computed by ZHETRF. */
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
/* >          as determined by ZHETRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N+NB+1,NB+3) */
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

/* > \date November 2015 */

/* > \ingroup complex16HEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhetri2x_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *nb, 
	integer *info, ftnlen uplo_len)
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
    extern /* Subroutine */ int zheswapr_(char *, integer *, doublecomplex *, 
	    integer *, integer *, integer *, ftnlen);
    static doublecomplex d__;
    static integer i__, j, k;
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


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 171 "zhetri2x.f"
    /* Parameter adjustments */
#line 171 "zhetri2x.f"
    a_dim1 = *lda;
#line 171 "zhetri2x.f"
    a_offset = 1 + a_dim1;
#line 171 "zhetri2x.f"
    a -= a_offset;
#line 171 "zhetri2x.f"
    --ipiv;
#line 171 "zhetri2x.f"
    work_dim1 = *n + *nb + 1;
#line 171 "zhetri2x.f"
    work_offset = 1 + work_dim1;
#line 171 "zhetri2x.f"
    work -= work_offset;
#line 171 "zhetri2x.f"

#line 171 "zhetri2x.f"
    /* Function Body */
#line 171 "zhetri2x.f"
    *info = 0;
#line 172 "zhetri2x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 173 "zhetri2x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 174 "zhetri2x.f"
	*info = -1;
#line 175 "zhetri2x.f"
    } else if (*n < 0) {
#line 176 "zhetri2x.f"
	*info = -2;
#line 177 "zhetri2x.f"
    } else if (*lda < max(1,*n)) {
#line 178 "zhetri2x.f"
	*info = -4;
#line 179 "zhetri2x.f"
    }

/*     Quick return if possible */


#line 184 "zhetri2x.f"
    if (*info != 0) {
#line 185 "zhetri2x.f"
	i__1 = -(*info);
#line 185 "zhetri2x.f"
	xerbla_("ZHETRI2X", &i__1, (ftnlen)8);
#line 186 "zhetri2x.f"
	return 0;
#line 187 "zhetri2x.f"
    }
#line 188 "zhetri2x.f"
    if (*n == 0) {
#line 188 "zhetri2x.f"
	return 0;
#line 188 "zhetri2x.f"
    }

/*     Convert A */
/*     Workspace got Non-diag elements of D */

#line 194 "zhetri2x.f"
    zsyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     Check that the diagonal matrix D is nonsingular. */

#line 198 "zhetri2x.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 202 "zhetri2x.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 203 "zhetri2x.f"
	    i__1 = *info + *info * a_dim1;
#line 203 "zhetri2x.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 203 "zhetri2x.f"
		return 0;
#line 203 "zhetri2x.f"
	    }
#line 205 "zhetri2x.f"
	}
#line 206 "zhetri2x.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 210 "zhetri2x.f"
	i__1 = *n;
#line 210 "zhetri2x.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 211 "zhetri2x.f"
	    i__2 = *info + *info * a_dim1;
#line 211 "zhetri2x.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 211 "zhetri2x.f"
		return 0;
#line 211 "zhetri2x.f"
	    }
#line 213 "zhetri2x.f"
	}
#line 214 "zhetri2x.f"
    }
#line 215 "zhetri2x.f"
    *info = 0;

/*  Splitting Workspace */
/*     U01 is a block (N,NB+1) */
/*     The first element of U01 is in WORK(1,1) */
/*     U11 is a block (NB+1,NB+1) */
/*     The first element of U11 is in WORK(N+1,1) */
#line 222 "zhetri2x.f"
    u11 = *n;
/*     INVD is a block (N,2) */
/*     The first element of INVD is in WORK(1,INVD) */
#line 225 "zhetri2x.f"
    invd = *nb + 2;
#line 227 "zhetri2x.f"
    if (upper) {

/*        invA = P * inv(U**H)*inv(D)*inv(U)*P**H. */

#line 231 "zhetri2x.f"
	ztrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 235 "zhetri2x.f"
	k = 1;
#line 236 "zhetri2x.f"
	while(k <= *n) {
#line 237 "zhetri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 239 "zhetri2x.f"
		i__1 = k + invd * work_dim1;
#line 239 "zhetri2x.f"
		i__2 = k + k * a_dim1;
#line 239 "zhetri2x.f"
		d__1 = 1. / a[i__2].r;
#line 239 "zhetri2x.f"
		work[i__1].r = d__1, work[i__1].i = 0.;
#line 240 "zhetri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 240 "zhetri2x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 241 "zhetri2x.f"
		++k;
#line 242 "zhetri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 244 "zhetri2x.f"
		d__1 = z_abs(&work[k + 1 + work_dim1]);
#line 244 "zhetri2x.f"
		t.r = d__1, t.i = 0.;
#line 245 "zhetri2x.f"
		i__1 = k + k * a_dim1;
#line 245 "zhetri2x.f"
		d__1 = a[i__1].r;
#line 245 "zhetri2x.f"
		z__2.r = d__1, z__2.i = 0.;
#line 245 "zhetri2x.f"
		z_div(&z__1, &z__2, &t);
#line 245 "zhetri2x.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 246 "zhetri2x.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 246 "zhetri2x.f"
		d__1 = a[i__1].r;
#line 246 "zhetri2x.f"
		z__2.r = d__1, z__2.i = 0.;
#line 246 "zhetri2x.f"
		z_div(&z__1, &z__2, &t);
#line 246 "zhetri2x.f"
		akp1.r = z__1.r, akp1.i = z__1.i;
#line 247 "zhetri2x.f"
		z_div(&z__1, &work[k + 1 + work_dim1], &t);
#line 247 "zhetri2x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 248 "zhetri2x.f"
		z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * 
			akp1.i + ak.i * akp1.r;
#line 248 "zhetri2x.f"
		z__2.r = z__3.r - 1., z__2.i = z__3.i;
#line 248 "zhetri2x.f"
		z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + 
			t.i * z__2.r;
#line 248 "zhetri2x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 249 "zhetri2x.f"
		i__1 = k + invd * work_dim1;
#line 249 "zhetri2x.f"
		z_div(&z__1, &akp1, &d__);
#line 249 "zhetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 250 "zhetri2x.f"
		i__1 = k + 1 + (invd + 1) * work_dim1;
#line 250 "zhetri2x.f"
		z_div(&z__1, &ak, &d__);
#line 250 "zhetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 251 "zhetri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 251 "zhetri2x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 251 "zhetri2x.f"
		z_div(&z__1, &z__2, &d__);
#line 251 "zhetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 252 "zhetri2x.f"
		i__1 = k + 1 + invd * work_dim1;
#line 252 "zhetri2x.f"
		d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
#line 252 "zhetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 253 "zhetri2x.f"
		k += 2;
#line 254 "zhetri2x.f"
	    }
#line 255 "zhetri2x.f"
	}

/*       inv(U**H) = (inv(U))**H */

/*       inv(U**H)*inv(D)*inv(U) */

#line 261 "zhetri2x.f"
	cut = *n;
#line 262 "zhetri2x.f"
	while(cut > 0) {
#line 263 "zhetri2x.f"
	    nnb = *nb;
#line 264 "zhetri2x.f"
	    if (cut <= nnb) {
#line 265 "zhetri2x.f"
		nnb = cut;
#line 266 "zhetri2x.f"
	    } else {
#line 267 "zhetri2x.f"
		count = 0;
/*             count negative elements, */
#line 269 "zhetri2x.f"
		i__1 = cut;
#line 269 "zhetri2x.f"
		for (i__ = cut + 1 - nnb; i__ <= i__1; ++i__) {
#line 270 "zhetri2x.f"
		    if (ipiv[i__] < 0) {
#line 270 "zhetri2x.f"
			++count;
#line 270 "zhetri2x.f"
		    }
#line 271 "zhetri2x.f"
		}
/*             need a even number for a clear cut */
#line 273 "zhetri2x.f"
		if (count % 2 == 1) {
#line 273 "zhetri2x.f"
		    ++nnb;
#line 273 "zhetri2x.f"
		}
#line 274 "zhetri2x.f"
	    }
#line 276 "zhetri2x.f"
	    cut -= nnb;

/*          U01 Block */

#line 280 "zhetri2x.f"
	    i__1 = cut;
#line 280 "zhetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 281 "zhetri2x.f"
		i__2 = nnb;
#line 281 "zhetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 282 "zhetri2x.f"
		    i__3 = i__ + j * work_dim1;
#line 282 "zhetri2x.f"
		    i__4 = i__ + (cut + j) * a_dim1;
#line 282 "zhetri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 283 "zhetri2x.f"
		}
#line 284 "zhetri2x.f"
	    }

/*          U11 Block */

#line 288 "zhetri2x.f"
	    i__1 = nnb;
#line 288 "zhetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 289 "zhetri2x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 289 "zhetri2x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 290 "zhetri2x.f"
		i__2 = i__ - 1;
#line 290 "zhetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 291 "zhetri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 291 "zhetri2x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 292 "zhetri2x.f"
		}
#line 293 "zhetri2x.f"
		i__2 = nnb;
#line 293 "zhetri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 294 "zhetri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 294 "zhetri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 294 "zhetri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 295 "zhetri2x.f"
		}
#line 296 "zhetri2x.f"
	    }

/*          invD*U01 */

#line 300 "zhetri2x.f"
	    i__ = 1;
#line 301 "zhetri2x.f"
	    while(i__ <= cut) {
#line 302 "zhetri2x.f"
		if (ipiv[i__] > 0) {
#line 303 "zhetri2x.f"
		    i__1 = nnb;
#line 303 "zhetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 304 "zhetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 304 "zhetri2x.f"
			i__3 = i__ + invd * work_dim1;
#line 304 "zhetri2x.f"
			i__4 = i__ + j * work_dim1;
#line 304 "zhetri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 304 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 305 "zhetri2x.f"
		    }
#line 306 "zhetri2x.f"
		    ++i__;
#line 307 "zhetri2x.f"
		} else {
#line 308 "zhetri2x.f"
		    i__1 = nnb;
#line 308 "zhetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 309 "zhetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 309 "zhetri2x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 310 "zhetri2x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 310 "zhetri2x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 311 "zhetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 311 "zhetri2x.f"
			i__3 = i__ + invd * work_dim1;
#line 311 "zhetri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 311 "zhetri2x.f"
			i__4 = i__ + (invd + 1) * work_dim1;
#line 311 "zhetri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 311 "zhetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 311 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 313 "zhetri2x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 313 "zhetri2x.f"
			i__3 = i__ + 1 + invd * work_dim1;
#line 313 "zhetri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 313 "zhetri2x.f"
			i__4 = i__ + 1 + (invd + 1) * work_dim1;
#line 313 "zhetri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 313 "zhetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 313 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 315 "zhetri2x.f"
		    }
#line 316 "zhetri2x.f"
		    i__ += 2;
#line 317 "zhetri2x.f"
		}
#line 318 "zhetri2x.f"
	    }

/*        invD1*U11 */

#line 322 "zhetri2x.f"
	    i__ = 1;
#line 323 "zhetri2x.f"
	    while(i__ <= nnb) {
#line 324 "zhetri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 325 "zhetri2x.f"
		    i__1 = nnb;
#line 325 "zhetri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 326 "zhetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 326 "zhetri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 326 "zhetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 326 "zhetri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 326 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 327 "zhetri2x.f"
		    }
#line 328 "zhetri2x.f"
		    ++i__;
#line 329 "zhetri2x.f"
		} else {
#line 330 "zhetri2x.f"
		    i__1 = nnb;
#line 330 "zhetri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 331 "zhetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 331 "zhetri2x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 332 "zhetri2x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 332 "zhetri2x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 333 "zhetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 333 "zhetri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 333 "zhetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 333 "zhetri2x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 333 "zhetri2x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 333 "zhetri2x.f"
			i__6 = u11 + i__ + 1 + j * work_dim1;
#line 333 "zhetri2x.f"
			z__3.r = work[i__5].r * work[i__6].r - work[i__5].i * 
				work[i__6].i, z__3.i = work[i__5].r * work[
				i__6].i + work[i__5].i * work[i__6].r;
#line 333 "zhetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 333 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 335 "zhetri2x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 335 "zhetri2x.f"
			i__3 = cut + i__ + 1 + invd * work_dim1;
#line 335 "zhetri2x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 335 "zhetri2x.f"
			i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
#line 335 "zhetri2x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 335 "zhetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 335 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 337 "zhetri2x.f"
		    }
#line 338 "zhetri2x.f"
		    i__ += 2;
#line 339 "zhetri2x.f"
		}
#line 340 "zhetri2x.f"
	    }

/*       U11**H*invD1*U11->U11 */

#line 344 "zhetri2x.f"
	    i__1 = *n + *nb + 1;
#line 344 "zhetri2x.f"
	    ztrmm_("L", "U", "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 
		    1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 347 "zhetri2x.f"
	    i__1 = nnb;
#line 347 "zhetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "zhetri2x.f"
		i__2 = nnb;
#line 348 "zhetri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 349 "zhetri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 349 "zhetri2x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 349 "zhetri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 350 "zhetri2x.f"
		}
#line 351 "zhetri2x.f"
	    }

/*          U01**H*invD*U01->A(CUT+I,CUT+J) */

#line 355 "zhetri2x.f"
	    i__1 = *n + *nb + 1;
#line 355 "zhetri2x.f"
	    i__2 = *n + *nb + 1;
#line 355 "zhetri2x.f"
	    zgemm_("C", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 
		    1], lda, &work[work_offset], &i__1, &c_b2, &work[u11 + 1 
		    + work_dim1], &i__2, (ftnlen)1, (ftnlen)1);

/*        U11 =  U11**H*invD1*U11 + U01**H*invD*U01 */

#line 360 "zhetri2x.f"
	    i__1 = nnb;
#line 360 "zhetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 361 "zhetri2x.f"
		i__2 = nnb;
#line 361 "zhetri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 362 "zhetri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 362 "zhetri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 362 "zhetri2x.f"
		    i__5 = u11 + i__ + j * work_dim1;
#line 362 "zhetri2x.f"
		    z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i + 
			    work[i__5].i;
#line 362 "zhetri2x.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 363 "zhetri2x.f"
		}
#line 364 "zhetri2x.f"
	    }

/*        U01 =  U00**H*invD0*U01 */

#line 368 "zhetri2x.f"
	    i__1 = *n + *nb + 1;
#line 368 "zhetri2x.f"
	    ztrmm_("L", uplo, "C", "U", &cut, &nnb, &c_b1, &a[a_offset], lda, 
		    &work[work_offset], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);

/*        Update U01 */

#line 374 "zhetri2x.f"
	    i__1 = cut;
#line 374 "zhetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 375 "zhetri2x.f"
		i__2 = nnb;
#line 375 "zhetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 376 "zhetri2x.f"
		    i__3 = i__ + (cut + j) * a_dim1;
#line 376 "zhetri2x.f"
		    i__4 = i__ + j * work_dim1;
#line 376 "zhetri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 377 "zhetri2x.f"
		}
#line 378 "zhetri2x.f"
	    }

/*      Next Block */

#line 382 "zhetri2x.f"
	}

/*        Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H */

#line 386 "zhetri2x.f"
	i__ = 1;
#line 387 "zhetri2x.f"
	while(i__ <= *n) {
#line 388 "zhetri2x.f"
	    if (ipiv[i__] > 0) {
#line 389 "zhetri2x.f"
		ip = ipiv[i__];
#line 390 "zhetri2x.f"
		if (i__ < ip) {
#line 390 "zhetri2x.f"
		    zheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 390 "zhetri2x.f"
		}
#line 391 "zhetri2x.f"
		if (i__ > ip) {
#line 391 "zhetri2x.f"
		    zheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 391 "zhetri2x.f"
		}
#line 392 "zhetri2x.f"
	    } else {
#line 393 "zhetri2x.f"
		ip = -ipiv[i__];
#line 394 "zhetri2x.f"
		++i__;
#line 395 "zhetri2x.f"
		if (i__ - 1 < ip) {
#line 395 "zhetri2x.f"
		    i__1 = i__ - 1;
#line 395 "zhetri2x.f"
		    zheswapr_(uplo, n, &a[a_offset], lda, &i__1, &ip, (ftnlen)
			    1);
#line 395 "zhetri2x.f"
		}
#line 397 "zhetri2x.f"
		if (i__ - 1 > ip) {
#line 397 "zhetri2x.f"
		    i__1 = i__ - 1;
#line 397 "zhetri2x.f"
		    zheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__1, (ftnlen)
			    1);
#line 397 "zhetri2x.f"
		}
#line 399 "zhetri2x.f"
	    }
#line 400 "zhetri2x.f"
	    ++i__;
#line 401 "zhetri2x.f"
	}
#line 402 "zhetri2x.f"
    } else {

/*        LOWER... */

/*        invA = P * inv(U**H)*inv(D)*inv(U)*P**H. */

#line 408 "zhetri2x.f"
	ztrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 412 "zhetri2x.f"
	k = *n;
#line 413 "zhetri2x.f"
	while(k >= 1) {
#line 414 "zhetri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 416 "zhetri2x.f"
		i__1 = k + invd * work_dim1;
#line 416 "zhetri2x.f"
		i__2 = k + k * a_dim1;
#line 416 "zhetri2x.f"
		d__1 = 1. / a[i__2].r;
#line 416 "zhetri2x.f"
		work[i__1].r = d__1, work[i__1].i = 0.;
#line 417 "zhetri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 417 "zhetri2x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 418 "zhetri2x.f"
		--k;
#line 419 "zhetri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 421 "zhetri2x.f"
		d__1 = z_abs(&work[k - 1 + work_dim1]);
#line 421 "zhetri2x.f"
		t.r = d__1, t.i = 0.;
#line 422 "zhetri2x.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 422 "zhetri2x.f"
		d__1 = a[i__1].r;
#line 422 "zhetri2x.f"
		z__2.r = d__1, z__2.i = 0.;
#line 422 "zhetri2x.f"
		z_div(&z__1, &z__2, &t);
#line 422 "zhetri2x.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 423 "zhetri2x.f"
		i__1 = k + k * a_dim1;
#line 423 "zhetri2x.f"
		d__1 = a[i__1].r;
#line 423 "zhetri2x.f"
		z__2.r = d__1, z__2.i = 0.;
#line 423 "zhetri2x.f"
		z_div(&z__1, &z__2, &t);
#line 423 "zhetri2x.f"
		akp1.r = z__1.r, akp1.i = z__1.i;
#line 424 "zhetri2x.f"
		z_div(&z__1, &work[k - 1 + work_dim1], &t);
#line 424 "zhetri2x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 425 "zhetri2x.f"
		z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * 
			akp1.i + ak.i * akp1.r;
#line 425 "zhetri2x.f"
		z__2.r = z__3.r - 1., z__2.i = z__3.i;
#line 425 "zhetri2x.f"
		z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + 
			t.i * z__2.r;
#line 425 "zhetri2x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 426 "zhetri2x.f"
		i__1 = k - 1 + invd * work_dim1;
#line 426 "zhetri2x.f"
		z_div(&z__1, &akp1, &d__);
#line 426 "zhetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 427 "zhetri2x.f"
		i__1 = k + invd * work_dim1;
#line 427 "zhetri2x.f"
		z_div(&z__1, &ak, &d__);
#line 427 "zhetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 428 "zhetri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 428 "zhetri2x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 428 "zhetri2x.f"
		z_div(&z__1, &z__2, &d__);
#line 428 "zhetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 429 "zhetri2x.f"
		i__1 = k - 1 + (invd + 1) * work_dim1;
#line 429 "zhetri2x.f"
		d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
#line 429 "zhetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 430 "zhetri2x.f"
		k += -2;
#line 431 "zhetri2x.f"
	    }
#line 432 "zhetri2x.f"
	}

/*       inv(U**H) = (inv(U))**H */

/*       inv(U**H)*inv(D)*inv(U) */

#line 438 "zhetri2x.f"
	cut = 0;
#line 439 "zhetri2x.f"
	while(cut < *n) {
#line 440 "zhetri2x.f"
	    nnb = *nb;
#line 441 "zhetri2x.f"
	    if (cut + nnb >= *n) {
#line 442 "zhetri2x.f"
		nnb = *n - cut;
#line 443 "zhetri2x.f"
	    } else {
#line 444 "zhetri2x.f"
		count = 0;
/*             count negative elements, */
#line 446 "zhetri2x.f"
		i__1 = cut + nnb;
#line 446 "zhetri2x.f"
		for (i__ = cut + 1; i__ <= i__1; ++i__) {
#line 447 "zhetri2x.f"
		    if (ipiv[i__] < 0) {
#line 447 "zhetri2x.f"
			++count;
#line 447 "zhetri2x.f"
		    }
#line 448 "zhetri2x.f"
		}
/*             need a even number for a clear cut */
#line 450 "zhetri2x.f"
		if (count % 2 == 1) {
#line 450 "zhetri2x.f"
		    ++nnb;
#line 450 "zhetri2x.f"
		}
#line 451 "zhetri2x.f"
	    }
/*      L21 Block */
#line 453 "zhetri2x.f"
	    i__1 = *n - cut - nnb;
#line 453 "zhetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 454 "zhetri2x.f"
		i__2 = nnb;
#line 454 "zhetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 455 "zhetri2x.f"
		    i__3 = i__ + j * work_dim1;
#line 455 "zhetri2x.f"
		    i__4 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 455 "zhetri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 456 "zhetri2x.f"
		}
#line 457 "zhetri2x.f"
	    }
/*     L11 Block */
#line 459 "zhetri2x.f"
	    i__1 = nnb;
#line 459 "zhetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 460 "zhetri2x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 460 "zhetri2x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 461 "zhetri2x.f"
		i__2 = nnb;
#line 461 "zhetri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 462 "zhetri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 462 "zhetri2x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 463 "zhetri2x.f"
		}
#line 464 "zhetri2x.f"
		i__2 = i__ - 1;
#line 464 "zhetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 465 "zhetri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 465 "zhetri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 465 "zhetri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 466 "zhetri2x.f"
		}
#line 467 "zhetri2x.f"
	    }

/*          invD*L21 */

#line 471 "zhetri2x.f"
	    i__ = *n - cut - nnb;
#line 472 "zhetri2x.f"
	    while(i__ >= 1) {
#line 473 "zhetri2x.f"
		if (ipiv[cut + nnb + i__] > 0) {
#line 474 "zhetri2x.f"
		    i__1 = nnb;
#line 474 "zhetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 475 "zhetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 475 "zhetri2x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 475 "zhetri2x.f"
			i__4 = i__ + j * work_dim1;
#line 475 "zhetri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 475 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 476 "zhetri2x.f"
		    }
#line 477 "zhetri2x.f"
		    --i__;
#line 478 "zhetri2x.f"
		} else {
#line 479 "zhetri2x.f"
		    i__1 = nnb;
#line 479 "zhetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 480 "zhetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 480 "zhetri2x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 481 "zhetri2x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 481 "zhetri2x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 482 "zhetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 482 "zhetri2x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 482 "zhetri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 482 "zhetri2x.f"
			i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
#line 482 "zhetri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 482 "zhetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 482 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 484 "zhetri2x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 484 "zhetri2x.f"
			i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
#line 484 "zhetri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 484 "zhetri2x.f"
			i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
#line 484 "zhetri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 484 "zhetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 484 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 486 "zhetri2x.f"
		    }
#line 487 "zhetri2x.f"
		    i__ += -2;
#line 488 "zhetri2x.f"
		}
#line 489 "zhetri2x.f"
	    }

/*        invD1*L11 */

#line 493 "zhetri2x.f"
	    i__ = nnb;
#line 494 "zhetri2x.f"
	    while(i__ >= 1) {
#line 495 "zhetri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 496 "zhetri2x.f"
		    i__1 = nnb;
#line 496 "zhetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 497 "zhetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 497 "zhetri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 497 "zhetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 497 "zhetri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 497 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 498 "zhetri2x.f"
		    }
#line 499 "zhetri2x.f"
		    --i__;
#line 500 "zhetri2x.f"
		} else {
#line 501 "zhetri2x.f"
		    i__1 = nnb;
#line 501 "zhetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 502 "zhetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 502 "zhetri2x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 503 "zhetri2x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 503 "zhetri2x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 504 "zhetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 504 "zhetri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 504 "zhetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 504 "zhetri2x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 504 "zhetri2x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 504 "zhetri2x.f"
			z__3.r = work[i__5].r * u11_ip1_j__.r - work[i__5].i *
				 u11_ip1_j__.i, z__3.i = work[i__5].r * 
				u11_ip1_j__.i + work[i__5].i * u11_ip1_j__.r;
#line 504 "zhetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 504 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 506 "zhetri2x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 506 "zhetri2x.f"
			i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
#line 506 "zhetri2x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 506 "zhetri2x.f"
			i__4 = cut + i__ - 1 + invd * work_dim1;
#line 506 "zhetri2x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 506 "zhetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 506 "zhetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 508 "zhetri2x.f"
		    }
#line 509 "zhetri2x.f"
		    i__ += -2;
#line 510 "zhetri2x.f"
		}
#line 511 "zhetri2x.f"
	    }

/*       L11**H*invD1*L11->L11 */

#line 515 "zhetri2x.f"
	    i__1 = *n + *nb + 1;
#line 515 "zhetri2x.f"
	    ztrmm_("L", uplo, "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 518 "zhetri2x.f"
	    i__1 = nnb;
#line 518 "zhetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 519 "zhetri2x.f"
		i__2 = i__;
#line 519 "zhetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 520 "zhetri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 520 "zhetri2x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 520 "zhetri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 521 "zhetri2x.f"
		}
#line 522 "zhetri2x.f"
	    }

#line 524 "zhetri2x.f"
	    if (cut + nnb < *n) {

/*          L21**H*invD2*L21->A(CUT+I,CUT+J) */

#line 528 "zhetri2x.f"
		i__1 = *n - nnb - cut;
#line 528 "zhetri2x.f"
		i__2 = *n + *nb + 1;
#line 528 "zhetri2x.f"
		i__3 = *n + *nb + 1;
#line 528 "zhetri2x.f"
		zgemm_("C", "N", &nnb, &nnb, &i__1, &c_b1, &a[cut + nnb + 1 + 
			(cut + 1) * a_dim1], lda, &work[work_offset], &i__2, &
			c_b2, &work[u11 + 1 + work_dim1], &i__3, (ftnlen)1, (
			ftnlen)1);

/*        L11 =  L11**H*invD1*L11 + U01**H*invD*U01 */

#line 534 "zhetri2x.f"
		i__1 = nnb;
#line 534 "zhetri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 535 "zhetri2x.f"
		    i__2 = i__;
#line 535 "zhetri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 536 "zhetri2x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 536 "zhetri2x.f"
			i__4 = cut + i__ + (cut + j) * a_dim1;
#line 536 "zhetri2x.f"
			i__5 = u11 + i__ + j * work_dim1;
#line 536 "zhetri2x.f"
			z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i 
				+ work[i__5].i;
#line 536 "zhetri2x.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 537 "zhetri2x.f"
		    }
#line 538 "zhetri2x.f"
		}

/*        L01 =  L22**H*invD2*L21 */

#line 542 "zhetri2x.f"
		i__1 = *n - nnb - cut;
#line 542 "zhetri2x.f"
		i__2 = *n + *nb + 1;
#line 542 "zhetri2x.f"
		ztrmm_("L", uplo, "C", "U", &i__1, &nnb, &c_b1, &a[cut + nnb 
			+ 1 + (cut + nnb + 1) * a_dim1], lda, &work[
			work_offset], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
/*      Update L21 */
#line 546 "zhetri2x.f"
		i__1 = *n - cut - nnb;
#line 546 "zhetri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 547 "zhetri2x.f"
		    i__2 = nnb;
#line 547 "zhetri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 548 "zhetri2x.f"
			i__3 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 548 "zhetri2x.f"
			i__4 = i__ + j * work_dim1;
#line 548 "zhetri2x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 549 "zhetri2x.f"
		    }
#line 550 "zhetri2x.f"
		}
#line 551 "zhetri2x.f"
	    } else {

/*        L11 =  L11**H*invD1*L11 */

#line 555 "zhetri2x.f"
		i__1 = nnb;
#line 555 "zhetri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 556 "zhetri2x.f"
		    i__2 = i__;
#line 556 "zhetri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 557 "zhetri2x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 557 "zhetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 557 "zhetri2x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 558 "zhetri2x.f"
		    }
#line 559 "zhetri2x.f"
		}
#line 560 "zhetri2x.f"
	    }

/*      Next Block */

#line 564 "zhetri2x.f"
	    cut += nnb;
#line 565 "zhetri2x.f"
	}

/*        Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H */

#line 569 "zhetri2x.f"
	i__ = *n;
#line 570 "zhetri2x.f"
	while(i__ >= 1) {
#line 571 "zhetri2x.f"
	    if (ipiv[i__] > 0) {
#line 572 "zhetri2x.f"
		ip = ipiv[i__];
#line 573 "zhetri2x.f"
		if (i__ < ip) {
#line 573 "zhetri2x.f"
		    zheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 573 "zhetri2x.f"
		}
#line 574 "zhetri2x.f"
		if (i__ > ip) {
#line 574 "zhetri2x.f"
		    zheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 574 "zhetri2x.f"
		}
#line 575 "zhetri2x.f"
	    } else {
#line 576 "zhetri2x.f"
		ip = -ipiv[i__];
#line 577 "zhetri2x.f"
		if (i__ < ip) {
#line 577 "zhetri2x.f"
		    zheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 577 "zhetri2x.f"
		}
#line 578 "zhetri2x.f"
		if (i__ > ip) {
#line 578 "zhetri2x.f"
		    zheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 578 "zhetri2x.f"
		}
#line 579 "zhetri2x.f"
		--i__;
#line 580 "zhetri2x.f"
	    }
#line 581 "zhetri2x.f"
	    --i__;
#line 582 "zhetri2x.f"
	}
#line 583 "zhetri2x.f"
    }

#line 585 "zhetri2x.f"
    return 0;

/*     End of ZHETRI2X */

} /* zhetri2x_ */

