#line 1 "chetri2x.f"
/* chetri2x.f -- translated by f2c (version 20100827).
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

#line 1 "chetri2x.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};

/* > \brief \b CHETRI2X */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRI2X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri2
x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri2
x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri2
x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), WORK( N+NB+1,* ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHETRI2X computes the inverse of a complex Hermitian indefinite matrix */
/* > A using the factorization A = U*D*U**H or A = L*D*L**H computed by */
/* > CHETRF. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the NNB diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by CHETRF. */
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
/* >          as determined by CHETRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N+NB+1,NB+3) */
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

/* > \ingroup complexHEcomputational */

/*  ===================================================================== */
/* Subroutine */ int chetri2x_(char *uplo, integer *n, doublecomplex *a, 
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
    extern /* Subroutine */ int cheswapr_(char *, integer *, doublecomplex *, 
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
    static integer iinfo;
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer count;
    static logical upper;
    static doublecomplex u01_i_j__, u11_i_j__;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), ctrtri_(
	    char *, char *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen, ftnlen), csyconv_(char *, char *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     ftnlen, ftnlen);
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

#line 171 "chetri2x.f"
    /* Parameter adjustments */
#line 171 "chetri2x.f"
    a_dim1 = *lda;
#line 171 "chetri2x.f"
    a_offset = 1 + a_dim1;
#line 171 "chetri2x.f"
    a -= a_offset;
#line 171 "chetri2x.f"
    --ipiv;
#line 171 "chetri2x.f"
    work_dim1 = *n + *nb + 1;
#line 171 "chetri2x.f"
    work_offset = 1 + work_dim1;
#line 171 "chetri2x.f"
    work -= work_offset;
#line 171 "chetri2x.f"

#line 171 "chetri2x.f"
    /* Function Body */
#line 171 "chetri2x.f"
    *info = 0;
#line 172 "chetri2x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 173 "chetri2x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 174 "chetri2x.f"
	*info = -1;
#line 175 "chetri2x.f"
    } else if (*n < 0) {
#line 176 "chetri2x.f"
	*info = -2;
#line 177 "chetri2x.f"
    } else if (*lda < max(1,*n)) {
#line 178 "chetri2x.f"
	*info = -4;
#line 179 "chetri2x.f"
    }

/*     Quick return if possible */


#line 184 "chetri2x.f"
    if (*info != 0) {
#line 185 "chetri2x.f"
	i__1 = -(*info);
#line 185 "chetri2x.f"
	xerbla_("CHETRI2X", &i__1, (ftnlen)8);
#line 186 "chetri2x.f"
	return 0;
#line 187 "chetri2x.f"
    }
#line 188 "chetri2x.f"
    if (*n == 0) {
#line 188 "chetri2x.f"
	return 0;
#line 188 "chetri2x.f"
    }

/*     Convert A */
/*     Workspace got Non-diag elements of D */

#line 194 "chetri2x.f"
    csyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     Check that the diagonal matrix D is nonsingular. */

#line 198 "chetri2x.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 202 "chetri2x.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 203 "chetri2x.f"
	    i__1 = *info + *info * a_dim1;
#line 203 "chetri2x.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 203 "chetri2x.f"
		return 0;
#line 203 "chetri2x.f"
	    }
#line 205 "chetri2x.f"
	}
#line 206 "chetri2x.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 210 "chetri2x.f"
	i__1 = *n;
#line 210 "chetri2x.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 211 "chetri2x.f"
	    i__2 = *info + *info * a_dim1;
#line 211 "chetri2x.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 211 "chetri2x.f"
		return 0;
#line 211 "chetri2x.f"
	    }
#line 213 "chetri2x.f"
	}
#line 214 "chetri2x.f"
    }
#line 215 "chetri2x.f"
    *info = 0;

/*  Splitting Workspace */
/*     U01 is a block (N,NB+1) */
/*     The first element of U01 is in WORK(1,1) */
/*     U11 is a block (NB+1,NB+1) */
/*     The first element of U11 is in WORK(N+1,1) */
#line 222 "chetri2x.f"
    u11 = *n;
/*     INVD is a block (N,2) */
/*     The first element of INVD is in WORK(1,INVD) */
#line 225 "chetri2x.f"
    invd = *nb + 2;
#line 227 "chetri2x.f"
    if (upper) {

/*        invA = P * inv(U**H)*inv(D)*inv(U)*P**H. */

#line 231 "chetri2x.f"
	ctrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 235 "chetri2x.f"
	k = 1;
#line 236 "chetri2x.f"
	while(k <= *n) {
#line 237 "chetri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 239 "chetri2x.f"
		i__1 = k + invd * work_dim1;
#line 239 "chetri2x.f"
		i__2 = k + k * a_dim1;
#line 239 "chetri2x.f"
		d__1 = 1. / a[i__2].r;
#line 239 "chetri2x.f"
		work[i__1].r = d__1, work[i__1].i = 0.;
#line 240 "chetri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 240 "chetri2x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 241 "chetri2x.f"
		++k;
#line 242 "chetri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 244 "chetri2x.f"
		d__1 = z_abs(&work[k + 1 + work_dim1]);
#line 244 "chetri2x.f"
		t.r = d__1, t.i = 0.;
#line 245 "chetri2x.f"
		i__1 = k + k * a_dim1;
#line 245 "chetri2x.f"
		d__1 = a[i__1].r;
#line 245 "chetri2x.f"
		z__2.r = d__1, z__2.i = 0.;
#line 245 "chetri2x.f"
		z_div(&z__1, &z__2, &t);
#line 245 "chetri2x.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 246 "chetri2x.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 246 "chetri2x.f"
		d__1 = a[i__1].r;
#line 246 "chetri2x.f"
		z__2.r = d__1, z__2.i = 0.;
#line 246 "chetri2x.f"
		z_div(&z__1, &z__2, &t);
#line 246 "chetri2x.f"
		akp1.r = z__1.r, akp1.i = z__1.i;
#line 247 "chetri2x.f"
		z_div(&z__1, &work[k + 1 + work_dim1], &t);
#line 247 "chetri2x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 248 "chetri2x.f"
		z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * 
			akp1.i + ak.i * akp1.r;
#line 248 "chetri2x.f"
		z__2.r = z__3.r - 1., z__2.i = z__3.i;
#line 248 "chetri2x.f"
		z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + 
			t.i * z__2.r;
#line 248 "chetri2x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 249 "chetri2x.f"
		i__1 = k + invd * work_dim1;
#line 249 "chetri2x.f"
		z_div(&z__1, &akp1, &d__);
#line 249 "chetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 250 "chetri2x.f"
		i__1 = k + 1 + (invd + 1) * work_dim1;
#line 250 "chetri2x.f"
		z_div(&z__1, &ak, &d__);
#line 250 "chetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 251 "chetri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 251 "chetri2x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 251 "chetri2x.f"
		z_div(&z__1, &z__2, &d__);
#line 251 "chetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 252 "chetri2x.f"
		i__1 = k + 1 + invd * work_dim1;
#line 252 "chetri2x.f"
		d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
#line 252 "chetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 253 "chetri2x.f"
		k += 2;
#line 254 "chetri2x.f"
	    }
#line 255 "chetri2x.f"
	}

/*       inv(U**H) = (inv(U))**H */

/*       inv(U**H)*inv(D)*inv(U) */

#line 261 "chetri2x.f"
	cut = *n;
#line 262 "chetri2x.f"
	while(cut > 0) {
#line 263 "chetri2x.f"
	    nnb = *nb;
#line 264 "chetri2x.f"
	    if (cut <= nnb) {
#line 265 "chetri2x.f"
		nnb = cut;
#line 266 "chetri2x.f"
	    } else {
#line 267 "chetri2x.f"
		count = 0;
/*             count negative elements, */
#line 269 "chetri2x.f"
		i__1 = cut;
#line 269 "chetri2x.f"
		for (i__ = cut + 1 - nnb; i__ <= i__1; ++i__) {
#line 270 "chetri2x.f"
		    if (ipiv[i__] < 0) {
#line 270 "chetri2x.f"
			++count;
#line 270 "chetri2x.f"
		    }
#line 271 "chetri2x.f"
		}
/*             need a even number for a clear cut */
#line 273 "chetri2x.f"
		if (count % 2 == 1) {
#line 273 "chetri2x.f"
		    ++nnb;
#line 273 "chetri2x.f"
		}
#line 274 "chetri2x.f"
	    }
#line 276 "chetri2x.f"
	    cut -= nnb;

/*          U01 Block */

#line 280 "chetri2x.f"
	    i__1 = cut;
#line 280 "chetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 281 "chetri2x.f"
		i__2 = nnb;
#line 281 "chetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 282 "chetri2x.f"
		    i__3 = i__ + j * work_dim1;
#line 282 "chetri2x.f"
		    i__4 = i__ + (cut + j) * a_dim1;
#line 282 "chetri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 283 "chetri2x.f"
		}
#line 284 "chetri2x.f"
	    }

/*          U11 Block */

#line 288 "chetri2x.f"
	    i__1 = nnb;
#line 288 "chetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 289 "chetri2x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 289 "chetri2x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 290 "chetri2x.f"
		i__2 = i__ - 1;
#line 290 "chetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 291 "chetri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 291 "chetri2x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 292 "chetri2x.f"
		}
#line 293 "chetri2x.f"
		i__2 = nnb;
#line 293 "chetri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 294 "chetri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 294 "chetri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 294 "chetri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 295 "chetri2x.f"
		}
#line 296 "chetri2x.f"
	    }

/*          invD*U01 */

#line 300 "chetri2x.f"
	    i__ = 1;
#line 301 "chetri2x.f"
	    while(i__ <= cut) {
#line 302 "chetri2x.f"
		if (ipiv[i__] > 0) {
#line 303 "chetri2x.f"
		    i__1 = nnb;
#line 303 "chetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 304 "chetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 304 "chetri2x.f"
			i__3 = i__ + invd * work_dim1;
#line 304 "chetri2x.f"
			i__4 = i__ + j * work_dim1;
#line 304 "chetri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 304 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 305 "chetri2x.f"
		    }
#line 306 "chetri2x.f"
		    ++i__;
#line 307 "chetri2x.f"
		} else {
#line 308 "chetri2x.f"
		    i__1 = nnb;
#line 308 "chetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 309 "chetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 309 "chetri2x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 310 "chetri2x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 310 "chetri2x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 311 "chetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 311 "chetri2x.f"
			i__3 = i__ + invd * work_dim1;
#line 311 "chetri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 311 "chetri2x.f"
			i__4 = i__ + (invd + 1) * work_dim1;
#line 311 "chetri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 311 "chetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 311 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 313 "chetri2x.f"
			i__2 = i__ + 1 + j * work_dim1;
#line 313 "chetri2x.f"
			i__3 = i__ + 1 + invd * work_dim1;
#line 313 "chetri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 313 "chetri2x.f"
			i__4 = i__ + 1 + (invd + 1) * work_dim1;
#line 313 "chetri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 313 "chetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 313 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 315 "chetri2x.f"
		    }
#line 316 "chetri2x.f"
		    i__ += 2;
#line 317 "chetri2x.f"
		}
#line 318 "chetri2x.f"
	    }

/*        invD1*U11 */

#line 322 "chetri2x.f"
	    i__ = 1;
#line 323 "chetri2x.f"
	    while(i__ <= nnb) {
#line 324 "chetri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 325 "chetri2x.f"
		    i__1 = nnb;
#line 325 "chetri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 326 "chetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 326 "chetri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 326 "chetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 326 "chetri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 326 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 327 "chetri2x.f"
		    }
#line 328 "chetri2x.f"
		    ++i__;
#line 329 "chetri2x.f"
		} else {
#line 330 "chetri2x.f"
		    i__1 = nnb;
#line 330 "chetri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 331 "chetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 331 "chetri2x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 332 "chetri2x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 332 "chetri2x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 333 "chetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 333 "chetri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 333 "chetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 333 "chetri2x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 333 "chetri2x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 333 "chetri2x.f"
			i__6 = u11 + i__ + 1 + j * work_dim1;
#line 333 "chetri2x.f"
			z__3.r = work[i__5].r * work[i__6].r - work[i__5].i * 
				work[i__6].i, z__3.i = work[i__5].r * work[
				i__6].i + work[i__5].i * work[i__6].r;
#line 333 "chetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 333 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 335 "chetri2x.f"
			i__2 = u11 + i__ + 1 + j * work_dim1;
#line 335 "chetri2x.f"
			i__3 = cut + i__ + 1 + invd * work_dim1;
#line 335 "chetri2x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 335 "chetri2x.f"
			i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
#line 335 "chetri2x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 335 "chetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 335 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 337 "chetri2x.f"
		    }
#line 338 "chetri2x.f"
		    i__ += 2;
#line 339 "chetri2x.f"
		}
#line 340 "chetri2x.f"
	    }

/*       U11**H*invD1*U11->U11 */

#line 344 "chetri2x.f"
	    i__1 = *n + *nb + 1;
#line 344 "chetri2x.f"
	    ctrmm_("L", "U", "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 
		    1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 347 "chetri2x.f"
	    i__1 = nnb;
#line 347 "chetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "chetri2x.f"
		i__2 = nnb;
#line 348 "chetri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 349 "chetri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 349 "chetri2x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 349 "chetri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 350 "chetri2x.f"
		}
#line 351 "chetri2x.f"
	    }

/*          U01**H*invD*U01->A(CUT+I,CUT+J) */

#line 355 "chetri2x.f"
	    i__1 = *n + *nb + 1;
#line 355 "chetri2x.f"
	    i__2 = *n + *nb + 1;
#line 355 "chetri2x.f"
	    cgemm_("C", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 
		    1], lda, &work[work_offset], &i__1, &c_b2, &work[u11 + 1 
		    + work_dim1], &i__2, (ftnlen)1, (ftnlen)1);

/*        U11 =  U11**H*invD1*U11 + U01**H*invD*U01 */

#line 360 "chetri2x.f"
	    i__1 = nnb;
#line 360 "chetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 361 "chetri2x.f"
		i__2 = nnb;
#line 361 "chetri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 362 "chetri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 362 "chetri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 362 "chetri2x.f"
		    i__5 = u11 + i__ + j * work_dim1;
#line 362 "chetri2x.f"
		    z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i + 
			    work[i__5].i;
#line 362 "chetri2x.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 363 "chetri2x.f"
		}
#line 364 "chetri2x.f"
	    }

/*        U01 =  U00**H*invD0*U01 */

#line 368 "chetri2x.f"
	    i__1 = *n + *nb + 1;
#line 368 "chetri2x.f"
	    ctrmm_("L", uplo, "C", "U", &cut, &nnb, &c_b1, &a[a_offset], lda, 
		    &work[work_offset], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);

/*        Update U01 */

#line 374 "chetri2x.f"
	    i__1 = cut;
#line 374 "chetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 375 "chetri2x.f"
		i__2 = nnb;
#line 375 "chetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 376 "chetri2x.f"
		    i__3 = i__ + (cut + j) * a_dim1;
#line 376 "chetri2x.f"
		    i__4 = i__ + j * work_dim1;
#line 376 "chetri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 377 "chetri2x.f"
		}
#line 378 "chetri2x.f"
	    }

/*      Next Block */

#line 382 "chetri2x.f"
	}

/*        Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H */

#line 386 "chetri2x.f"
	i__ = 1;
#line 387 "chetri2x.f"
	while(i__ <= *n) {
#line 388 "chetri2x.f"
	    if (ipiv[i__] > 0) {
#line 389 "chetri2x.f"
		ip = ipiv[i__];
#line 390 "chetri2x.f"
		if (i__ < ip) {
#line 390 "chetri2x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 390 "chetri2x.f"
		}
#line 391 "chetri2x.f"
		if (i__ > ip) {
#line 391 "chetri2x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 391 "chetri2x.f"
		}
#line 392 "chetri2x.f"
	    } else {
#line 393 "chetri2x.f"
		ip = -ipiv[i__];
#line 394 "chetri2x.f"
		++i__;
#line 395 "chetri2x.f"
		if (i__ - 1 < ip) {
#line 395 "chetri2x.f"
		    i__1 = i__ - 1;
#line 395 "chetri2x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &i__1, &ip, (ftnlen)
			    1);
#line 395 "chetri2x.f"
		}
#line 397 "chetri2x.f"
		if (i__ - 1 > ip) {
#line 397 "chetri2x.f"
		    i__1 = i__ - 1;
#line 397 "chetri2x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__1, (ftnlen)
			    1);
#line 397 "chetri2x.f"
		}
#line 399 "chetri2x.f"
	    }
#line 400 "chetri2x.f"
	    ++i__;
#line 401 "chetri2x.f"
	}
#line 402 "chetri2x.f"
    } else {

/*        LOWER... */

/*        invA = P * inv(U**H)*inv(D)*inv(U)*P**H. */

#line 408 "chetri2x.f"
	ctrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 412 "chetri2x.f"
	k = *n;
#line 413 "chetri2x.f"
	while(k >= 1) {
#line 414 "chetri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 416 "chetri2x.f"
		i__1 = k + invd * work_dim1;
#line 416 "chetri2x.f"
		i__2 = k + k * a_dim1;
#line 416 "chetri2x.f"
		d__1 = 1. / a[i__2].r;
#line 416 "chetri2x.f"
		work[i__1].r = d__1, work[i__1].i = 0.;
#line 417 "chetri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 417 "chetri2x.f"
		work[i__1].r = 0., work[i__1].i = 0.;
#line 418 "chetri2x.f"
		--k;
#line 419 "chetri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 421 "chetri2x.f"
		d__1 = z_abs(&work[k - 1 + work_dim1]);
#line 421 "chetri2x.f"
		t.r = d__1, t.i = 0.;
#line 422 "chetri2x.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 422 "chetri2x.f"
		d__1 = a[i__1].r;
#line 422 "chetri2x.f"
		z__2.r = d__1, z__2.i = 0.;
#line 422 "chetri2x.f"
		z_div(&z__1, &z__2, &t);
#line 422 "chetri2x.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 423 "chetri2x.f"
		i__1 = k + k * a_dim1;
#line 423 "chetri2x.f"
		d__1 = a[i__1].r;
#line 423 "chetri2x.f"
		z__2.r = d__1, z__2.i = 0.;
#line 423 "chetri2x.f"
		z_div(&z__1, &z__2, &t);
#line 423 "chetri2x.f"
		akp1.r = z__1.r, akp1.i = z__1.i;
#line 424 "chetri2x.f"
		z_div(&z__1, &work[k - 1 + work_dim1], &t);
#line 424 "chetri2x.f"
		akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 425 "chetri2x.f"
		z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * 
			akp1.i + ak.i * akp1.r;
#line 425 "chetri2x.f"
		z__2.r = z__3.r - 1., z__2.i = z__3.i;
#line 425 "chetri2x.f"
		z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + 
			t.i * z__2.r;
#line 425 "chetri2x.f"
		d__.r = z__1.r, d__.i = z__1.i;
#line 426 "chetri2x.f"
		i__1 = k - 1 + invd * work_dim1;
#line 426 "chetri2x.f"
		z_div(&z__1, &akp1, &d__);
#line 426 "chetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 427 "chetri2x.f"
		i__1 = k + invd * work_dim1;
#line 427 "chetri2x.f"
		z_div(&z__1, &ak, &d__);
#line 427 "chetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 428 "chetri2x.f"
		i__1 = k + (invd + 1) * work_dim1;
#line 428 "chetri2x.f"
		z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 428 "chetri2x.f"
		z_div(&z__1, &z__2, &d__);
#line 428 "chetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 429 "chetri2x.f"
		i__1 = k - 1 + (invd + 1) * work_dim1;
#line 429 "chetri2x.f"
		d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
#line 429 "chetri2x.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 430 "chetri2x.f"
		k += -2;
#line 431 "chetri2x.f"
	    }
#line 432 "chetri2x.f"
	}

/*       inv(U**H) = (inv(U))**H */

/*       inv(U**H)*inv(D)*inv(U) */

#line 438 "chetri2x.f"
	cut = 0;
#line 439 "chetri2x.f"
	while(cut < *n) {
#line 440 "chetri2x.f"
	    nnb = *nb;
#line 441 "chetri2x.f"
	    if (cut + nnb >= *n) {
#line 442 "chetri2x.f"
		nnb = *n - cut;
#line 443 "chetri2x.f"
	    } else {
#line 444 "chetri2x.f"
		count = 0;
/*             count negative elements, */
#line 446 "chetri2x.f"
		i__1 = cut + nnb;
#line 446 "chetri2x.f"
		for (i__ = cut + 1; i__ <= i__1; ++i__) {
#line 447 "chetri2x.f"
		    if (ipiv[i__] < 0) {
#line 447 "chetri2x.f"
			++count;
#line 447 "chetri2x.f"
		    }
#line 448 "chetri2x.f"
		}
/*             need a even number for a clear cut */
#line 450 "chetri2x.f"
		if (count % 2 == 1) {
#line 450 "chetri2x.f"
		    ++nnb;
#line 450 "chetri2x.f"
		}
#line 451 "chetri2x.f"
	    }
/*      L21 Block */
#line 453 "chetri2x.f"
	    i__1 = *n - cut - nnb;
#line 453 "chetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 454 "chetri2x.f"
		i__2 = nnb;
#line 454 "chetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 455 "chetri2x.f"
		    i__3 = i__ + j * work_dim1;
#line 455 "chetri2x.f"
		    i__4 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 455 "chetri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 456 "chetri2x.f"
		}
#line 457 "chetri2x.f"
	    }
/*     L11 Block */
#line 459 "chetri2x.f"
	    i__1 = nnb;
#line 459 "chetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 460 "chetri2x.f"
		i__2 = u11 + i__ + i__ * work_dim1;
#line 460 "chetri2x.f"
		work[i__2].r = 1., work[i__2].i = 0.;
#line 461 "chetri2x.f"
		i__2 = nnb;
#line 461 "chetri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 462 "chetri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 462 "chetri2x.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 463 "chetri2x.f"
		}
#line 464 "chetri2x.f"
		i__2 = i__ - 1;
#line 464 "chetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 465 "chetri2x.f"
		    i__3 = u11 + i__ + j * work_dim1;
#line 465 "chetri2x.f"
		    i__4 = cut + i__ + (cut + j) * a_dim1;
#line 465 "chetri2x.f"
		    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
#line 466 "chetri2x.f"
		}
#line 467 "chetri2x.f"
	    }

/*          invD*L21 */

#line 471 "chetri2x.f"
	    i__ = *n - cut - nnb;
#line 472 "chetri2x.f"
	    while(i__ >= 1) {
#line 473 "chetri2x.f"
		if (ipiv[cut + nnb + i__] > 0) {
#line 474 "chetri2x.f"
		    i__1 = nnb;
#line 474 "chetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 475 "chetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 475 "chetri2x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 475 "chetri2x.f"
			i__4 = i__ + j * work_dim1;
#line 475 "chetri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 475 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 476 "chetri2x.f"
		    }
#line 477 "chetri2x.f"
		    --i__;
#line 478 "chetri2x.f"
		} else {
#line 479 "chetri2x.f"
		    i__1 = nnb;
#line 479 "chetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 480 "chetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 480 "chetri2x.f"
			u01_i_j__.r = work[i__2].r, u01_i_j__.i = work[i__2]
				.i;
#line 481 "chetri2x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 481 "chetri2x.f"
			u01_ip1_j__.r = work[i__2].r, u01_ip1_j__.i = work[
				i__2].i;
#line 482 "chetri2x.f"
			i__2 = i__ + j * work_dim1;
#line 482 "chetri2x.f"
			i__3 = cut + nnb + i__ + invd * work_dim1;
#line 482 "chetri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 482 "chetri2x.f"
			i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
#line 482 "chetri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 482 "chetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 482 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 484 "chetri2x.f"
			i__2 = i__ - 1 + j * work_dim1;
#line 484 "chetri2x.f"
			i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
#line 484 "chetri2x.f"
			z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * 
				u01_i_j__.i, z__2.i = work[i__3].r * 
				u01_i_j__.i + work[i__3].i * u01_i_j__.r;
#line 484 "chetri2x.f"
			i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
#line 484 "chetri2x.f"
			z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i *
				 u01_ip1_j__.i, z__3.i = work[i__4].r * 
				u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r;
#line 484 "chetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 484 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 486 "chetri2x.f"
		    }
#line 487 "chetri2x.f"
		    i__ += -2;
#line 488 "chetri2x.f"
		}
#line 489 "chetri2x.f"
	    }

/*        invD1*L11 */

#line 493 "chetri2x.f"
	    i__ = nnb;
#line 494 "chetri2x.f"
	    while(i__ >= 1) {
#line 495 "chetri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 496 "chetri2x.f"
		    i__1 = nnb;
#line 496 "chetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 497 "chetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 497 "chetri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 497 "chetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 497 "chetri2x.f"
			z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__1.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 497 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 498 "chetri2x.f"
		    }
#line 499 "chetri2x.f"
		    --i__;
#line 500 "chetri2x.f"
		} else {
#line 501 "chetri2x.f"
		    i__1 = nnb;
#line 501 "chetri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 502 "chetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 502 "chetri2x.f"
			u11_i_j__.r = work[i__2].r, u11_i_j__.i = work[i__2]
				.i;
#line 503 "chetri2x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 503 "chetri2x.f"
			u11_ip1_j__.r = work[i__2].r, u11_ip1_j__.i = work[
				i__2].i;
#line 504 "chetri2x.f"
			i__2 = u11 + i__ + j * work_dim1;
#line 504 "chetri2x.f"
			i__3 = cut + i__ + invd * work_dim1;
#line 504 "chetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 504 "chetri2x.f"
			z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * 
				work[i__4].i, z__2.i = work[i__3].r * work[
				i__4].i + work[i__3].i * work[i__4].r;
#line 504 "chetri2x.f"
			i__5 = cut + i__ + (invd + 1) * work_dim1;
#line 504 "chetri2x.f"
			z__3.r = work[i__5].r * u11_ip1_j__.r - work[i__5].i *
				 u11_ip1_j__.i, z__3.i = work[i__5].r * 
				u11_ip1_j__.i + work[i__5].i * u11_ip1_j__.r;
#line 504 "chetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 504 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 506 "chetri2x.f"
			i__2 = u11 + i__ - 1 + j * work_dim1;
#line 506 "chetri2x.f"
			i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
#line 506 "chetri2x.f"
			z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * 
				u11_i_j__.i, z__2.i = work[i__3].r * 
				u11_i_j__.i + work[i__3].i * u11_i_j__.r;
#line 506 "chetri2x.f"
			i__4 = cut + i__ - 1 + invd * work_dim1;
#line 506 "chetri2x.f"
			z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i *
				 u11_ip1_j__.i, z__3.i = work[i__4].r * 
				u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r;
#line 506 "chetri2x.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 506 "chetri2x.f"
			work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 508 "chetri2x.f"
		    }
#line 509 "chetri2x.f"
		    i__ += -2;
#line 510 "chetri2x.f"
		}
#line 511 "chetri2x.f"
	    }

/*       L11**H*invD1*L11->L11 */

#line 515 "chetri2x.f"
	    i__1 = *n + *nb + 1;
#line 515 "chetri2x.f"
	    ctrmm_("L", uplo, "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 518 "chetri2x.f"
	    i__1 = nnb;
#line 518 "chetri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 519 "chetri2x.f"
		i__2 = i__;
#line 519 "chetri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 520 "chetri2x.f"
		    i__3 = cut + i__ + (cut + j) * a_dim1;
#line 520 "chetri2x.f"
		    i__4 = u11 + i__ + j * work_dim1;
#line 520 "chetri2x.f"
		    a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 521 "chetri2x.f"
		}
#line 522 "chetri2x.f"
	    }

#line 524 "chetri2x.f"
	    if (cut + nnb < *n) {

/*          L21**H*invD2*L21->A(CUT+I,CUT+J) */

#line 528 "chetri2x.f"
		i__1 = *n - nnb - cut;
#line 528 "chetri2x.f"
		i__2 = *n + *nb + 1;
#line 528 "chetri2x.f"
		i__3 = *n + *nb + 1;
#line 528 "chetri2x.f"
		cgemm_("C", "N", &nnb, &nnb, &i__1, &c_b1, &a[cut + nnb + 1 + 
			(cut + 1) * a_dim1], lda, &work[work_offset], &i__2, &
			c_b2, &work[u11 + 1 + work_dim1], &i__3, (ftnlen)1, (
			ftnlen)1);

/*        L11 =  L11**H*invD1*L11 + U01**H*invD*U01 */

#line 534 "chetri2x.f"
		i__1 = nnb;
#line 534 "chetri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 535 "chetri2x.f"
		    i__2 = i__;
#line 535 "chetri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 536 "chetri2x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 536 "chetri2x.f"
			i__4 = cut + i__ + (cut + j) * a_dim1;
#line 536 "chetri2x.f"
			i__5 = u11 + i__ + j * work_dim1;
#line 536 "chetri2x.f"
			z__1.r = a[i__4].r + work[i__5].r, z__1.i = a[i__4].i 
				+ work[i__5].i;
#line 536 "chetri2x.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 537 "chetri2x.f"
		    }
#line 538 "chetri2x.f"
		}

/*        L01 =  L22**H*invD2*L21 */

#line 542 "chetri2x.f"
		i__1 = *n - nnb - cut;
#line 542 "chetri2x.f"
		i__2 = *n + *nb + 1;
#line 542 "chetri2x.f"
		ctrmm_("L", uplo, "C", "U", &i__1, &nnb, &c_b1, &a[cut + nnb 
			+ 1 + (cut + nnb + 1) * a_dim1], lda, &work[
			work_offset], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
/*      Update L21 */
#line 546 "chetri2x.f"
		i__1 = *n - cut - nnb;
#line 546 "chetri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 547 "chetri2x.f"
		    i__2 = nnb;
#line 547 "chetri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 548 "chetri2x.f"
			i__3 = cut + nnb + i__ + (cut + j) * a_dim1;
#line 548 "chetri2x.f"
			i__4 = i__ + j * work_dim1;
#line 548 "chetri2x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 549 "chetri2x.f"
		    }
#line 550 "chetri2x.f"
		}
#line 551 "chetri2x.f"
	    } else {

/*        L11 =  L11**H*invD1*L11 */

#line 555 "chetri2x.f"
		i__1 = nnb;
#line 555 "chetri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 556 "chetri2x.f"
		    i__2 = i__;
#line 556 "chetri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 557 "chetri2x.f"
			i__3 = cut + i__ + (cut + j) * a_dim1;
#line 557 "chetri2x.f"
			i__4 = u11 + i__ + j * work_dim1;
#line 557 "chetri2x.f"
			a[i__3].r = work[i__4].r, a[i__3].i = work[i__4].i;
#line 558 "chetri2x.f"
		    }
#line 559 "chetri2x.f"
		}
#line 560 "chetri2x.f"
	    }

/*      Next Block */

#line 564 "chetri2x.f"
	    cut += nnb;
#line 565 "chetri2x.f"
	}

/*        Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H */

#line 569 "chetri2x.f"
	i__ = *n;
#line 570 "chetri2x.f"
	while(i__ >= 1) {
#line 571 "chetri2x.f"
	    if (ipiv[i__] > 0) {
#line 572 "chetri2x.f"
		ip = ipiv[i__];
#line 573 "chetri2x.f"
		if (i__ < ip) {
#line 573 "chetri2x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 573 "chetri2x.f"
		}
#line 574 "chetri2x.f"
		if (i__ > ip) {
#line 574 "chetri2x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 574 "chetri2x.f"
		}
#line 575 "chetri2x.f"
	    } else {
#line 576 "chetri2x.f"
		ip = -ipiv[i__];
#line 577 "chetri2x.f"
		if (i__ < ip) {
#line 577 "chetri2x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 577 "chetri2x.f"
		}
#line 578 "chetri2x.f"
		if (i__ > ip) {
#line 578 "chetri2x.f"
		    cheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 578 "chetri2x.f"
		}
#line 579 "chetri2x.f"
		--i__;
#line 580 "chetri2x.f"
	    }
#line 581 "chetri2x.f"
	    --i__;
#line 582 "chetri2x.f"
	}
#line 583 "chetri2x.f"
    }

#line 585 "chetri2x.f"
    return 0;

/*     End of CHETRI2X */

} /* chetri2x_ */

