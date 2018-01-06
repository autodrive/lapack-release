#line 1 "ssytri2x.f"
/* ssytri2x.f -- translated by f2c (version 20100827).
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

#line 1 "ssytri2x.f"
/* Table of constant values */

static doublereal c_b11 = 1.;
static doublereal c_b15 = 0.;

/* > \brief \b SSYTRI2X */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTRI2X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri2
x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri2
x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri2
x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), WORK( N+NB+1,* ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTRI2X computes the inverse of a real symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > SSYTRF. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the NNB diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by SSYTRF. */
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
/* >          as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N+NNB+1,NNB+3) */
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

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssytri2x_(char *uplo, integer *n, doublereal *a, integer 
	*lda, integer *ipiv, doublereal *work, integer *nb, integer *info, 
	ftnlen uplo_len)
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
    static integer iinfo;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer count;
    static logical upper;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal u01_i_j__, u11_i_j__;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), strtri_(
	    char *, char *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static doublereal u01_ip1_j__, u11_ip1_j__;
    extern /* Subroutine */ int ssyconv_(char *, char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
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

#line 168 "ssytri2x.f"
    /* Parameter adjustments */
#line 168 "ssytri2x.f"
    a_dim1 = *lda;
#line 168 "ssytri2x.f"
    a_offset = 1 + a_dim1;
#line 168 "ssytri2x.f"
    a -= a_offset;
#line 168 "ssytri2x.f"
    --ipiv;
#line 168 "ssytri2x.f"
    work_dim1 = *n + *nb + 1;
#line 168 "ssytri2x.f"
    work_offset = 1 + work_dim1;
#line 168 "ssytri2x.f"
    work -= work_offset;
#line 168 "ssytri2x.f"

#line 168 "ssytri2x.f"
    /* Function Body */
#line 168 "ssytri2x.f"
    *info = 0;
#line 169 "ssytri2x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 170 "ssytri2x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 171 "ssytri2x.f"
	*info = -1;
#line 172 "ssytri2x.f"
    } else if (*n < 0) {
#line 173 "ssytri2x.f"
	*info = -2;
#line 174 "ssytri2x.f"
    } else if (*lda < max(1,*n)) {
#line 175 "ssytri2x.f"
	*info = -4;
#line 176 "ssytri2x.f"
    }

/*     Quick return if possible */


#line 181 "ssytri2x.f"
    if (*info != 0) {
#line 182 "ssytri2x.f"
	i__1 = -(*info);
#line 182 "ssytri2x.f"
	xerbla_("SSYTRI2X", &i__1, (ftnlen)8);
#line 183 "ssytri2x.f"
	return 0;
#line 184 "ssytri2x.f"
    }
#line 185 "ssytri2x.f"
    if (*n == 0) {
#line 185 "ssytri2x.f"
	return 0;
#line 185 "ssytri2x.f"
    }

/*     Convert A */
/*     Workspace got Non-diag elements of D */

#line 191 "ssytri2x.f"
    ssyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     Check that the diagonal matrix D is nonsingular. */

#line 195 "ssytri2x.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 199 "ssytri2x.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 200 "ssytri2x.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 200 "ssytri2x.f"
		return 0;
#line 200 "ssytri2x.f"
	    }
#line 202 "ssytri2x.f"
	}
#line 203 "ssytri2x.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 207 "ssytri2x.f"
	i__1 = *n;
#line 207 "ssytri2x.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 208 "ssytri2x.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 208 "ssytri2x.f"
		return 0;
#line 208 "ssytri2x.f"
	    }
#line 210 "ssytri2x.f"
	}
#line 211 "ssytri2x.f"
    }
#line 212 "ssytri2x.f"
    *info = 0;

/*  Splitting Workspace */
/*     U01 is a block (N,NB+1) */
/*     The first element of U01 is in WORK(1,1) */
/*     U11 is a block (NB+1,NB+1) */
/*     The first element of U11 is in WORK(N+1,1) */
#line 219 "ssytri2x.f"
    u11 = *n;
/*     INVD is a block (N,2) */
/*     The first element of INVD is in WORK(1,INVD) */
#line 222 "ssytri2x.f"
    invd = *nb + 2;
#line 224 "ssytri2x.f"
    if (upper) {

/*        invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */

#line 228 "ssytri2x.f"
	strtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 232 "ssytri2x.f"
	k = 1;
#line 233 "ssytri2x.f"
	while(k <= *n) {
#line 234 "ssytri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 236 "ssytri2x.f"
		work[k + invd * work_dim1] = 1. / a[k + k * a_dim1];
#line 237 "ssytri2x.f"
		work[k + (invd + 1) * work_dim1] = 0.;
#line 238 "ssytri2x.f"
		++k;
#line 239 "ssytri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 241 "ssytri2x.f"
		t = work[k + 1 + work_dim1];
#line 242 "ssytri2x.f"
		ak = a[k + k * a_dim1] / t;
#line 243 "ssytri2x.f"
		akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
#line 244 "ssytri2x.f"
		akkp1 = work[k + 1 + work_dim1] / t;
#line 245 "ssytri2x.f"
		d__ = t * (ak * akp1 - 1.);
#line 246 "ssytri2x.f"
		work[k + invd * work_dim1] = akp1 / d__;
#line 247 "ssytri2x.f"
		work[k + 1 + (invd + 1) * work_dim1] = ak / d__;
#line 248 "ssytri2x.f"
		work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
#line 249 "ssytri2x.f"
		work[k + 1 + invd * work_dim1] = -akkp1 / d__;
#line 250 "ssytri2x.f"
		k += 2;
#line 251 "ssytri2x.f"
	    }
#line 252 "ssytri2x.f"
	}

/*       inv(U**T) = (inv(U))**T */

/*       inv(U**T)*inv(D)*inv(U) */

#line 258 "ssytri2x.f"
	cut = *n;
#line 259 "ssytri2x.f"
	while(cut > 0) {
#line 260 "ssytri2x.f"
	    nnb = *nb;
#line 261 "ssytri2x.f"
	    if (cut <= nnb) {
#line 262 "ssytri2x.f"
		nnb = cut;
#line 263 "ssytri2x.f"
	    } else {
#line 264 "ssytri2x.f"
		count = 0;
/*             count negative elements, */
#line 266 "ssytri2x.f"
		i__1 = cut;
#line 266 "ssytri2x.f"
		for (i__ = cut + 1 - nnb; i__ <= i__1; ++i__) {
#line 267 "ssytri2x.f"
		    if (ipiv[i__] < 0) {
#line 267 "ssytri2x.f"
			++count;
#line 267 "ssytri2x.f"
		    }
#line 268 "ssytri2x.f"
		}
/*             need a even number for a clear cut */
#line 270 "ssytri2x.f"
		if (count % 2 == 1) {
#line 270 "ssytri2x.f"
		    ++nnb;
#line 270 "ssytri2x.f"
		}
#line 271 "ssytri2x.f"
	    }
#line 273 "ssytri2x.f"
	    cut -= nnb;

/*          U01 Block */

#line 277 "ssytri2x.f"
	    i__1 = cut;
#line 277 "ssytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 278 "ssytri2x.f"
		i__2 = nnb;
#line 278 "ssytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 279 "ssytri2x.f"
		    work[i__ + j * work_dim1] = a[i__ + (cut + j) * a_dim1];
#line 280 "ssytri2x.f"
		}
#line 281 "ssytri2x.f"
	    }

/*          U11 Block */

#line 285 "ssytri2x.f"
	    i__1 = nnb;
#line 285 "ssytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 286 "ssytri2x.f"
		work[u11 + i__ + i__ * work_dim1] = 1.;
#line 287 "ssytri2x.f"
		i__2 = i__ - 1;
#line 287 "ssytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 288 "ssytri2x.f"
		    work[u11 + i__ + j * work_dim1] = 0.;
#line 289 "ssytri2x.f"
		}
#line 290 "ssytri2x.f"
		i__2 = nnb;
#line 290 "ssytri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 291 "ssytri2x.f"
		    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) 
			    * a_dim1];
#line 292 "ssytri2x.f"
		}
#line 293 "ssytri2x.f"
	    }

/*          invD*U01 */

#line 297 "ssytri2x.f"
	    i__ = 1;
#line 298 "ssytri2x.f"
	    while(i__ <= cut) {
#line 299 "ssytri2x.f"
		if (ipiv[i__] > 0) {
#line 300 "ssytri2x.f"
		    i__1 = nnb;
#line 300 "ssytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 301 "ssytri2x.f"
			work[i__ + j * work_dim1] = work[i__ + invd * 
				work_dim1] * work[i__ + j * work_dim1];
#line 302 "ssytri2x.f"
		    }
#line 303 "ssytri2x.f"
		    ++i__;
#line 304 "ssytri2x.f"
		} else {
#line 305 "ssytri2x.f"
		    i__1 = nnb;
#line 305 "ssytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 306 "ssytri2x.f"
			u01_i_j__ = work[i__ + j * work_dim1];
#line 307 "ssytri2x.f"
			u01_ip1_j__ = work[i__ + 1 + j * work_dim1];
#line 308 "ssytri2x.f"
			work[i__ + j * work_dim1] = work[i__ + invd * 
				work_dim1] * u01_i_j__ + work[i__ + (invd + 1)
				 * work_dim1] * u01_ip1_j__;
#line 310 "ssytri2x.f"
			work[i__ + 1 + j * work_dim1] = work[i__ + 1 + invd * 
				work_dim1] * u01_i_j__ + work[i__ + 1 + (invd 
				+ 1) * work_dim1] * u01_ip1_j__;
#line 312 "ssytri2x.f"
		    }
#line 313 "ssytri2x.f"
		    i__ += 2;
#line 314 "ssytri2x.f"
		}
#line 315 "ssytri2x.f"
	    }

/*        invD1*U11 */

#line 319 "ssytri2x.f"
	    i__ = 1;
#line 320 "ssytri2x.f"
	    while(i__ <= nnb) {
#line 321 "ssytri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 322 "ssytri2x.f"
		    i__1 = nnb;
#line 322 "ssytri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 323 "ssytri2x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1];
#line 324 "ssytri2x.f"
		    }
#line 325 "ssytri2x.f"
		    ++i__;
#line 326 "ssytri2x.f"
		} else {
#line 327 "ssytri2x.f"
		    i__1 = nnb;
#line 327 "ssytri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 328 "ssytri2x.f"
			u11_i_j__ = work[u11 + i__ + j * work_dim1];
#line 329 "ssytri2x.f"
			u11_ip1_j__ = work[u11 + i__ + 1 + j * work_dim1];
#line 330 "ssytri2x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1] + work[cut + i__ + (invd + 1) * 
				work_dim1] * work[u11 + i__ + 1 + j * 
				work_dim1];
#line 332 "ssytri2x.f"
			work[u11 + i__ + 1 + j * work_dim1] = work[cut + i__ 
				+ 1 + invd * work_dim1] * u11_i_j__ + work[
				cut + i__ + 1 + (invd + 1) * work_dim1] * 
				u11_ip1_j__;
#line 334 "ssytri2x.f"
		    }
#line 335 "ssytri2x.f"
		    i__ += 2;
#line 336 "ssytri2x.f"
		}
#line 337 "ssytri2x.f"
	    }

/*       U11**T*invD1*U11->U11 */

#line 341 "ssytri2x.f"
	    i__1 = *n + *nb + 1;
#line 341 "ssytri2x.f"
	    strmm_("L", "U", "T", "U", &nnb, &nnb, &c_b11, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 344 "ssytri2x.f"
	    i__1 = nnb;
#line 344 "ssytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 345 "ssytri2x.f"
		i__2 = nnb;
#line 345 "ssytri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 346 "ssytri2x.f"
		    a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * 
			    work_dim1];
#line 347 "ssytri2x.f"
		}
#line 348 "ssytri2x.f"
	    }

/*          U01**T*invD*U01->A(CUT+I,CUT+J) */

#line 352 "ssytri2x.f"
	    i__1 = *n + *nb + 1;
#line 352 "ssytri2x.f"
	    i__2 = *n + *nb + 1;
#line 352 "ssytri2x.f"
	    sgemm_("T", "N", &nnb, &nnb, &cut, &c_b11, &a[(cut + 1) * a_dim1 
		    + 1], lda, &work[work_offset], &i__1, &c_b15, &work[u11 + 
		    1 + work_dim1], &i__2, (ftnlen)1, (ftnlen)1);

/*        U11 =  U11**T*invD1*U11 + U01**T*invD*U01 */

#line 357 "ssytri2x.f"
	    i__1 = nnb;
#line 357 "ssytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 358 "ssytri2x.f"
		i__2 = nnb;
#line 358 "ssytri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 359 "ssytri2x.f"
		    a[cut + i__ + (cut + j) * a_dim1] += work[u11 + i__ + j * 
			    work_dim1];
#line 360 "ssytri2x.f"
		}
#line 361 "ssytri2x.f"
	    }

/*        U01 =  U00**T*invD0*U01 */

#line 365 "ssytri2x.f"
	    i__1 = *n + *nb + 1;
#line 365 "ssytri2x.f"
	    strmm_("L", uplo, "T", "U", &cut, &nnb, &c_b11, &a[a_offset], lda,
		     &work[work_offset], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);

/*        Update U01 */

#line 371 "ssytri2x.f"
	    i__1 = cut;
#line 371 "ssytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 372 "ssytri2x.f"
		i__2 = nnb;
#line 372 "ssytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 373 "ssytri2x.f"
		    a[i__ + (cut + j) * a_dim1] = work[i__ + j * work_dim1];
#line 374 "ssytri2x.f"
		}
#line 375 "ssytri2x.f"
	    }

/*      Next Block */

#line 379 "ssytri2x.f"
	}

/*        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */

#line 383 "ssytri2x.f"
	i__ = 1;
#line 384 "ssytri2x.f"
	while(i__ <= *n) {
#line 385 "ssytri2x.f"
	    if (ipiv[i__] > 0) {
#line 386 "ssytri2x.f"
		ip = ipiv[i__];
#line 387 "ssytri2x.f"
		if (i__ < ip) {
#line 387 "ssytri2x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 387 "ssytri2x.f"
		}
#line 388 "ssytri2x.f"
		if (i__ > ip) {
#line 388 "ssytri2x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 388 "ssytri2x.f"
		}
#line 389 "ssytri2x.f"
	    } else {
#line 390 "ssytri2x.f"
		ip = -ipiv[i__];
#line 391 "ssytri2x.f"
		++i__;
#line 392 "ssytri2x.f"
		if (i__ - 1 < ip) {
#line 392 "ssytri2x.f"
		    i__1 = i__ - 1;
#line 392 "ssytri2x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &i__1, &ip, (ftnlen)
			    1);
#line 392 "ssytri2x.f"
		}
#line 394 "ssytri2x.f"
		if (i__ - 1 > ip) {
#line 394 "ssytri2x.f"
		    i__1 = i__ - 1;
#line 394 "ssytri2x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__1, (ftnlen)
			    1);
#line 394 "ssytri2x.f"
		}
#line 396 "ssytri2x.f"
	    }
#line 397 "ssytri2x.f"
	    ++i__;
#line 398 "ssytri2x.f"
	}
#line 399 "ssytri2x.f"
    } else {

/*        LOWER... */

/*        invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */

#line 405 "ssytri2x.f"
	strtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 409 "ssytri2x.f"
	k = *n;
#line 410 "ssytri2x.f"
	while(k >= 1) {
#line 411 "ssytri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 413 "ssytri2x.f"
		work[k + invd * work_dim1] = 1. / a[k + k * a_dim1];
#line 414 "ssytri2x.f"
		work[k + (invd + 1) * work_dim1] = 0.;
#line 415 "ssytri2x.f"
		--k;
#line 416 "ssytri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 418 "ssytri2x.f"
		t = work[k - 1 + work_dim1];
#line 419 "ssytri2x.f"
		ak = a[k - 1 + (k - 1) * a_dim1] / t;
#line 420 "ssytri2x.f"
		akp1 = a[k + k * a_dim1] / t;
#line 421 "ssytri2x.f"
		akkp1 = work[k - 1 + work_dim1] / t;
#line 422 "ssytri2x.f"
		d__ = t * (ak * akp1 - 1.);
#line 423 "ssytri2x.f"
		work[k - 1 + invd * work_dim1] = akp1 / d__;
#line 424 "ssytri2x.f"
		work[k + invd * work_dim1] = ak / d__;
#line 425 "ssytri2x.f"
		work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
#line 426 "ssytri2x.f"
		work[k - 1 + (invd + 1) * work_dim1] = -akkp1 / d__;
#line 427 "ssytri2x.f"
		k += -2;
#line 428 "ssytri2x.f"
	    }
#line 429 "ssytri2x.f"
	}

/*       inv(U**T) = (inv(U))**T */

/*       inv(U**T)*inv(D)*inv(U) */

#line 435 "ssytri2x.f"
	cut = 0;
#line 436 "ssytri2x.f"
	while(cut < *n) {
#line 437 "ssytri2x.f"
	    nnb = *nb;
#line 438 "ssytri2x.f"
	    if (cut + nnb > *n) {
#line 439 "ssytri2x.f"
		nnb = *n - cut;
#line 440 "ssytri2x.f"
	    } else {
#line 441 "ssytri2x.f"
		count = 0;
/*             count negative elements, */
#line 443 "ssytri2x.f"
		i__1 = cut + nnb;
#line 443 "ssytri2x.f"
		for (i__ = cut + 1; i__ <= i__1; ++i__) {
#line 444 "ssytri2x.f"
		    if (ipiv[i__] < 0) {
#line 444 "ssytri2x.f"
			++count;
#line 444 "ssytri2x.f"
		    }
#line 445 "ssytri2x.f"
		}
/*             need a even number for a clear cut */
#line 447 "ssytri2x.f"
		if (count % 2 == 1) {
#line 447 "ssytri2x.f"
		    ++nnb;
#line 447 "ssytri2x.f"
		}
#line 448 "ssytri2x.f"
	    }
/*     L21 Block */
#line 450 "ssytri2x.f"
	    i__1 = *n - cut - nnb;
#line 450 "ssytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 451 "ssytri2x.f"
		i__2 = nnb;
#line 451 "ssytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 452 "ssytri2x.f"
		    work[i__ + j * work_dim1] = a[cut + nnb + i__ + (cut + j) 
			    * a_dim1];
#line 453 "ssytri2x.f"
		}
#line 454 "ssytri2x.f"
	    }
/*     L11 Block */
#line 456 "ssytri2x.f"
	    i__1 = nnb;
#line 456 "ssytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 457 "ssytri2x.f"
		work[u11 + i__ + i__ * work_dim1] = 1.;
#line 458 "ssytri2x.f"
		i__2 = nnb;
#line 458 "ssytri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 459 "ssytri2x.f"
		    work[u11 + i__ + j * work_dim1] = 0.;
#line 460 "ssytri2x.f"
		}
#line 461 "ssytri2x.f"
		i__2 = i__ - 1;
#line 461 "ssytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 462 "ssytri2x.f"
		    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) 
			    * a_dim1];
#line 463 "ssytri2x.f"
		}
#line 464 "ssytri2x.f"
	    }

/*          invD*L21 */

#line 468 "ssytri2x.f"
	    i__ = *n - cut - nnb;
#line 469 "ssytri2x.f"
	    while(i__ >= 1) {
#line 470 "ssytri2x.f"
		if (ipiv[cut + nnb + i__] > 0) {
#line 471 "ssytri2x.f"
		    i__1 = nnb;
#line 471 "ssytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 472 "ssytri2x.f"
			work[i__ + j * work_dim1] = work[cut + nnb + i__ + 
				invd * work_dim1] * work[i__ + j * work_dim1];
#line 473 "ssytri2x.f"
		    }
#line 474 "ssytri2x.f"
		    --i__;
#line 475 "ssytri2x.f"
		} else {
#line 476 "ssytri2x.f"
		    i__1 = nnb;
#line 476 "ssytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 477 "ssytri2x.f"
			u01_i_j__ = work[i__ + j * work_dim1];
#line 478 "ssytri2x.f"
			u01_ip1_j__ = work[i__ - 1 + j * work_dim1];
#line 479 "ssytri2x.f"
			work[i__ + j * work_dim1] = work[cut + nnb + i__ + 
				invd * work_dim1] * u01_i_j__ + work[cut + 
				nnb + i__ + (invd + 1) * work_dim1] * 
				u01_ip1_j__;
#line 481 "ssytri2x.f"
			work[i__ - 1 + j * work_dim1] = work[cut + nnb + i__ 
				- 1 + (invd + 1) * work_dim1] * u01_i_j__ + 
				work[cut + nnb + i__ - 1 + invd * work_dim1] *
				 u01_ip1_j__;
#line 483 "ssytri2x.f"
		    }
#line 484 "ssytri2x.f"
		    i__ += -2;
#line 485 "ssytri2x.f"
		}
#line 486 "ssytri2x.f"
	    }

/*        invD1*L11 */

#line 490 "ssytri2x.f"
	    i__ = nnb;
#line 491 "ssytri2x.f"
	    while(i__ >= 1) {
#line 492 "ssytri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 493 "ssytri2x.f"
		    i__1 = nnb;
#line 493 "ssytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 494 "ssytri2x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1];
#line 495 "ssytri2x.f"
		    }
#line 496 "ssytri2x.f"
		    --i__;
#line 497 "ssytri2x.f"
		} else {
#line 498 "ssytri2x.f"
		    i__1 = nnb;
#line 498 "ssytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 499 "ssytri2x.f"
			u11_i_j__ = work[u11 + i__ + j * work_dim1];
#line 500 "ssytri2x.f"
			u11_ip1_j__ = work[u11 + i__ - 1 + j * work_dim1];
#line 501 "ssytri2x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1] + work[cut + i__ + (invd + 1) * 
				work_dim1] * u11_ip1_j__;
#line 503 "ssytri2x.f"
			work[u11 + i__ - 1 + j * work_dim1] = work[cut + i__ 
				- 1 + (invd + 1) * work_dim1] * u11_i_j__ + 
				work[cut + i__ - 1 + invd * work_dim1] * 
				u11_ip1_j__;
#line 505 "ssytri2x.f"
		    }
#line 506 "ssytri2x.f"
		    i__ += -2;
#line 507 "ssytri2x.f"
		}
#line 508 "ssytri2x.f"
	    }

/*       L11**T*invD1*L11->L11 */

#line 512 "ssytri2x.f"
	    i__1 = *n + *nb + 1;
#line 512 "ssytri2x.f"
	    strmm_("L", uplo, "T", "U", &nnb, &nnb, &c_b11, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 516 "ssytri2x.f"
	    i__1 = nnb;
#line 516 "ssytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 517 "ssytri2x.f"
		i__2 = i__;
#line 517 "ssytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 518 "ssytri2x.f"
		    a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * 
			    work_dim1];
#line 519 "ssytri2x.f"
		}
#line 520 "ssytri2x.f"
	    }

#line 522 "ssytri2x.f"
	    if (cut + nnb < *n) {

/*          L21**T*invD2*L21->A(CUT+I,CUT+J) */

#line 526 "ssytri2x.f"
		i__1 = *n - nnb - cut;
#line 526 "ssytri2x.f"
		i__2 = *n + *nb + 1;
#line 526 "ssytri2x.f"
		i__3 = *n + *nb + 1;
#line 526 "ssytri2x.f"
		sgemm_("T", "N", &nnb, &nnb, &i__1, &c_b11, &a[cut + nnb + 1 
			+ (cut + 1) * a_dim1], lda, &work[work_offset], &i__2,
			 &c_b15, &work[u11 + 1 + work_dim1], &i__3, (ftnlen)1,
			 (ftnlen)1);

/*        L11 =  L11**T*invD1*L11 + U01**T*invD*U01 */

#line 532 "ssytri2x.f"
		i__1 = nnb;
#line 532 "ssytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 533 "ssytri2x.f"
		    i__2 = i__;
#line 533 "ssytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 534 "ssytri2x.f"
			a[cut + i__ + (cut + j) * a_dim1] += work[u11 + i__ + 
				j * work_dim1];
#line 535 "ssytri2x.f"
		    }
#line 536 "ssytri2x.f"
		}

/*        L01 =  L22**T*invD2*L21 */

#line 540 "ssytri2x.f"
		i__1 = *n - nnb - cut;
#line 540 "ssytri2x.f"
		i__2 = *n + *nb + 1;
#line 540 "ssytri2x.f"
		strmm_("L", uplo, "T", "U", &i__1, &nnb, &c_b11, &a[cut + nnb 
			+ 1 + (cut + nnb + 1) * a_dim1], lda, &work[
			work_offset], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);

/*      Update L21 */

#line 545 "ssytri2x.f"
		i__1 = *n - cut - nnb;
#line 545 "ssytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 546 "ssytri2x.f"
		    i__2 = nnb;
#line 546 "ssytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 547 "ssytri2x.f"
			a[cut + nnb + i__ + (cut + j) * a_dim1] = work[i__ + 
				j * work_dim1];
#line 548 "ssytri2x.f"
		    }
#line 549 "ssytri2x.f"
		}
#line 551 "ssytri2x.f"
	    } else {

/*        L11 =  L11**T*invD1*L11 */

#line 555 "ssytri2x.f"
		i__1 = nnb;
#line 555 "ssytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 556 "ssytri2x.f"
		    i__2 = i__;
#line 556 "ssytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 557 "ssytri2x.f"
			a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + 
				j * work_dim1];
#line 558 "ssytri2x.f"
		    }
#line 559 "ssytri2x.f"
		}
#line 560 "ssytri2x.f"
	    }

/*      Next Block */

#line 564 "ssytri2x.f"
	    cut += nnb;
#line 565 "ssytri2x.f"
	}

/*        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */

#line 569 "ssytri2x.f"
	i__ = *n;
#line 570 "ssytri2x.f"
	while(i__ >= 1) {
#line 571 "ssytri2x.f"
	    if (ipiv[i__] > 0) {
#line 572 "ssytri2x.f"
		ip = ipiv[i__];
#line 573 "ssytri2x.f"
		if (i__ < ip) {
#line 573 "ssytri2x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 573 "ssytri2x.f"
		}
#line 574 "ssytri2x.f"
		if (i__ > ip) {
#line 574 "ssytri2x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 574 "ssytri2x.f"
		}
#line 575 "ssytri2x.f"
	    } else {
#line 576 "ssytri2x.f"
		ip = -ipiv[i__];
#line 577 "ssytri2x.f"
		if (i__ < ip) {
#line 577 "ssytri2x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 577 "ssytri2x.f"
		}
#line 578 "ssytri2x.f"
		if (i__ > ip) {
#line 578 "ssytri2x.f"
		    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 578 "ssytri2x.f"
		}
#line 579 "ssytri2x.f"
		--i__;
#line 580 "ssytri2x.f"
	    }
#line 581 "ssytri2x.f"
	    --i__;
#line 582 "ssytri2x.f"
	}
#line 583 "ssytri2x.f"
    }

#line 585 "ssytri2x.f"
    return 0;

/*     End of SSYTRI2X */

} /* ssytri2x_ */

