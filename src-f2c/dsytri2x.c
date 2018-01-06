#line 1 "dsytri2x.f"
/* dsytri2x.f -- translated by f2c (version 20100827).
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

#line 1 "dsytri2x.f"
/* Table of constant values */

static doublereal c_b11 = 1.;
static doublereal c_b15 = 0.;

/* > \brief \b DSYTRI2X */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRI2X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytri2
x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytri2
x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytri2
x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N, NB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( N+NB+1,* ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTRI2X computes the inverse of a real symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > DSYTRF. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the NNB diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by DSYTRF. */
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
/* >          as determined by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N+NNB+1,NNB+3) */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsytri2x_(char *uplo, integer *n, doublereal *a, integer 
	*lda, integer *ipiv, doublereal *work, integer *nb, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int dsyswapr_(char *, integer *, doublereal *, 
	    integer *, integer *, integer *, ftnlen);
    static doublereal d__;
    static integer i__, j, k;
    static doublereal t, ak;
    static integer u11, ip, nnb, cut;
    static doublereal akp1;
    static integer invd;
    static doublereal akkp1;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer count;
    static logical upper;
    static doublereal u01_i_j__, u11_i_j__;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dtrtri_(
	    char *, char *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dsyconv_(char *, char *, integer *, doublereal *,
	     integer *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal u01_ip1_j__, u11_ip1_j__;


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

#line 168 "dsytri2x.f"
    /* Parameter adjustments */
#line 168 "dsytri2x.f"
    a_dim1 = *lda;
#line 168 "dsytri2x.f"
    a_offset = 1 + a_dim1;
#line 168 "dsytri2x.f"
    a -= a_offset;
#line 168 "dsytri2x.f"
    --ipiv;
#line 168 "dsytri2x.f"
    work_dim1 = *n + *nb + 1;
#line 168 "dsytri2x.f"
    work_offset = 1 + work_dim1;
#line 168 "dsytri2x.f"
    work -= work_offset;
#line 168 "dsytri2x.f"

#line 168 "dsytri2x.f"
    /* Function Body */
#line 168 "dsytri2x.f"
    *info = 0;
#line 169 "dsytri2x.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 170 "dsytri2x.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 171 "dsytri2x.f"
	*info = -1;
#line 172 "dsytri2x.f"
    } else if (*n < 0) {
#line 173 "dsytri2x.f"
	*info = -2;
#line 174 "dsytri2x.f"
    } else if (*lda < max(1,*n)) {
#line 175 "dsytri2x.f"
	*info = -4;
#line 176 "dsytri2x.f"
    }

/*     Quick return if possible */


#line 181 "dsytri2x.f"
    if (*info != 0) {
#line 182 "dsytri2x.f"
	i__1 = -(*info);
#line 182 "dsytri2x.f"
	xerbla_("DSYTRI2X", &i__1, (ftnlen)8);
#line 183 "dsytri2x.f"
	return 0;
#line 184 "dsytri2x.f"
    }
#line 185 "dsytri2x.f"
    if (*n == 0) {
#line 185 "dsytri2x.f"
	return 0;
#line 185 "dsytri2x.f"
    }

/*     Convert A */
/*     Workspace got Non-diag elements of D */

#line 191 "dsytri2x.f"
    dsyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], &
	    iinfo, (ftnlen)1, (ftnlen)1);

/*     Check that the diagonal matrix D is nonsingular. */

#line 195 "dsytri2x.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 199 "dsytri2x.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 200 "dsytri2x.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 200 "dsytri2x.f"
		return 0;
#line 200 "dsytri2x.f"
	    }
#line 202 "dsytri2x.f"
	}
#line 203 "dsytri2x.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 207 "dsytri2x.f"
	i__1 = *n;
#line 207 "dsytri2x.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 208 "dsytri2x.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 208 "dsytri2x.f"
		return 0;
#line 208 "dsytri2x.f"
	    }
#line 210 "dsytri2x.f"
	}
#line 211 "dsytri2x.f"
    }
#line 212 "dsytri2x.f"
    *info = 0;

/*  Splitting Workspace */
/*     U01 is a block (N,NB+1) */
/*     The first element of U01 is in WORK(1,1) */
/*     U11 is a block (NB+1,NB+1) */
/*     The first element of U11 is in WORK(N+1,1) */
#line 219 "dsytri2x.f"
    u11 = *n;
/*     INVD is a block (N,2) */
/*     The first element of INVD is in WORK(1,INVD) */
#line 222 "dsytri2x.f"
    invd = *nb + 2;
#line 224 "dsytri2x.f"
    if (upper) {

/*        invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */

#line 228 "dsytri2x.f"
	dtrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 232 "dsytri2x.f"
	k = 1;
#line 233 "dsytri2x.f"
	while(k <= *n) {
#line 234 "dsytri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 236 "dsytri2x.f"
		work[k + invd * work_dim1] = 1. / a[k + k * a_dim1];
#line 237 "dsytri2x.f"
		work[k + (invd + 1) * work_dim1] = 0.;
#line 238 "dsytri2x.f"
		++k;
#line 239 "dsytri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 241 "dsytri2x.f"
		t = work[k + 1 + work_dim1];
#line 242 "dsytri2x.f"
		ak = a[k + k * a_dim1] / t;
#line 243 "dsytri2x.f"
		akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
#line 244 "dsytri2x.f"
		akkp1 = work[k + 1 + work_dim1] / t;
#line 245 "dsytri2x.f"
		d__ = t * (ak * akp1 - 1.);
#line 246 "dsytri2x.f"
		work[k + invd * work_dim1] = akp1 / d__;
#line 247 "dsytri2x.f"
		work[k + 1 + (invd + 1) * work_dim1] = ak / d__;
#line 248 "dsytri2x.f"
		work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
#line 249 "dsytri2x.f"
		work[k + 1 + invd * work_dim1] = -akkp1 / d__;
#line 250 "dsytri2x.f"
		k += 2;
#line 251 "dsytri2x.f"
	    }
#line 252 "dsytri2x.f"
	}

/*       inv(U**T) = (inv(U))**T */

/*       inv(U**T)*inv(D)*inv(U) */

#line 258 "dsytri2x.f"
	cut = *n;
#line 259 "dsytri2x.f"
	while(cut > 0) {
#line 260 "dsytri2x.f"
	    nnb = *nb;
#line 261 "dsytri2x.f"
	    if (cut <= nnb) {
#line 262 "dsytri2x.f"
		nnb = cut;
#line 263 "dsytri2x.f"
	    } else {
#line 264 "dsytri2x.f"
		count = 0;
/*             count negative elements, */
#line 266 "dsytri2x.f"
		i__1 = cut;
#line 266 "dsytri2x.f"
		for (i__ = cut + 1 - nnb; i__ <= i__1; ++i__) {
#line 267 "dsytri2x.f"
		    if (ipiv[i__] < 0) {
#line 267 "dsytri2x.f"
			++count;
#line 267 "dsytri2x.f"
		    }
#line 268 "dsytri2x.f"
		}
/*             need a even number for a clear cut */
#line 270 "dsytri2x.f"
		if (count % 2 == 1) {
#line 270 "dsytri2x.f"
		    ++nnb;
#line 270 "dsytri2x.f"
		}
#line 271 "dsytri2x.f"
	    }
#line 273 "dsytri2x.f"
	    cut -= nnb;

/*          U01 Block */

#line 277 "dsytri2x.f"
	    i__1 = cut;
#line 277 "dsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 278 "dsytri2x.f"
		i__2 = nnb;
#line 278 "dsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 279 "dsytri2x.f"
		    work[i__ + j * work_dim1] = a[i__ + (cut + j) * a_dim1];
#line 280 "dsytri2x.f"
		}
#line 281 "dsytri2x.f"
	    }

/*          U11 Block */

#line 285 "dsytri2x.f"
	    i__1 = nnb;
#line 285 "dsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 286 "dsytri2x.f"
		work[u11 + i__ + i__ * work_dim1] = 1.;
#line 287 "dsytri2x.f"
		i__2 = i__ - 1;
#line 287 "dsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 288 "dsytri2x.f"
		    work[u11 + i__ + j * work_dim1] = 0.;
#line 289 "dsytri2x.f"
		}
#line 290 "dsytri2x.f"
		i__2 = nnb;
#line 290 "dsytri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 291 "dsytri2x.f"
		    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) 
			    * a_dim1];
#line 292 "dsytri2x.f"
		}
#line 293 "dsytri2x.f"
	    }

/*          invD*U01 */

#line 297 "dsytri2x.f"
	    i__ = 1;
#line 298 "dsytri2x.f"
	    while(i__ <= cut) {
#line 299 "dsytri2x.f"
		if (ipiv[i__] > 0) {
#line 300 "dsytri2x.f"
		    i__1 = nnb;
#line 300 "dsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 301 "dsytri2x.f"
			work[i__ + j * work_dim1] = work[i__ + invd * 
				work_dim1] * work[i__ + j * work_dim1];
#line 302 "dsytri2x.f"
		    }
#line 303 "dsytri2x.f"
		    ++i__;
#line 304 "dsytri2x.f"
		} else {
#line 305 "dsytri2x.f"
		    i__1 = nnb;
#line 305 "dsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 306 "dsytri2x.f"
			u01_i_j__ = work[i__ + j * work_dim1];
#line 307 "dsytri2x.f"
			u01_ip1_j__ = work[i__ + 1 + j * work_dim1];
#line 308 "dsytri2x.f"
			work[i__ + j * work_dim1] = work[i__ + invd * 
				work_dim1] * u01_i_j__ + work[i__ + (invd + 1)
				 * work_dim1] * u01_ip1_j__;
#line 310 "dsytri2x.f"
			work[i__ + 1 + j * work_dim1] = work[i__ + 1 + invd * 
				work_dim1] * u01_i_j__ + work[i__ + 1 + (invd 
				+ 1) * work_dim1] * u01_ip1_j__;
#line 312 "dsytri2x.f"
		    }
#line 313 "dsytri2x.f"
		    i__ += 2;
#line 314 "dsytri2x.f"
		}
#line 315 "dsytri2x.f"
	    }

/*        invD1*U11 */

#line 319 "dsytri2x.f"
	    i__ = 1;
#line 320 "dsytri2x.f"
	    while(i__ <= nnb) {
#line 321 "dsytri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 322 "dsytri2x.f"
		    i__1 = nnb;
#line 322 "dsytri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 323 "dsytri2x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1];
#line 324 "dsytri2x.f"
		    }
#line 325 "dsytri2x.f"
		    ++i__;
#line 326 "dsytri2x.f"
		} else {
#line 327 "dsytri2x.f"
		    i__1 = nnb;
#line 327 "dsytri2x.f"
		    for (j = i__; j <= i__1; ++j) {
#line 328 "dsytri2x.f"
			u11_i_j__ = work[u11 + i__ + j * work_dim1];
#line 329 "dsytri2x.f"
			u11_ip1_j__ = work[u11 + i__ + 1 + j * work_dim1];
#line 330 "dsytri2x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1] + work[cut + i__ + (invd + 1) * 
				work_dim1] * work[u11 + i__ + 1 + j * 
				work_dim1];
#line 332 "dsytri2x.f"
			work[u11 + i__ + 1 + j * work_dim1] = work[cut + i__ 
				+ 1 + invd * work_dim1] * u11_i_j__ + work[
				cut + i__ + 1 + (invd + 1) * work_dim1] * 
				u11_ip1_j__;
#line 334 "dsytri2x.f"
		    }
#line 335 "dsytri2x.f"
		    i__ += 2;
#line 336 "dsytri2x.f"
		}
#line 337 "dsytri2x.f"
	    }

/*       U11**T*invD1*U11->U11 */

#line 341 "dsytri2x.f"
	    i__1 = *n + *nb + 1;
#line 341 "dsytri2x.f"
	    dtrmm_("L", "U", "T", "U", &nnb, &nnb, &c_b11, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 344 "dsytri2x.f"
	    i__1 = nnb;
#line 344 "dsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 345 "dsytri2x.f"
		i__2 = nnb;
#line 345 "dsytri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 346 "dsytri2x.f"
		    a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * 
			    work_dim1];
#line 347 "dsytri2x.f"
		}
#line 348 "dsytri2x.f"
	    }

/*          U01**T*invD*U01->A(CUT+I,CUT+J) */

#line 352 "dsytri2x.f"
	    i__1 = *n + *nb + 1;
#line 352 "dsytri2x.f"
	    i__2 = *n + *nb + 1;
#line 352 "dsytri2x.f"
	    dgemm_("T", "N", &nnb, &nnb, &cut, &c_b11, &a[(cut + 1) * a_dim1 
		    + 1], lda, &work[work_offset], &i__1, &c_b15, &work[u11 + 
		    1 + work_dim1], &i__2, (ftnlen)1, (ftnlen)1);

/*        U11 =  U11**T*invD1*U11 + U01**T*invD*U01 */

#line 358 "dsytri2x.f"
	    i__1 = nnb;
#line 358 "dsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 359 "dsytri2x.f"
		i__2 = nnb;
#line 359 "dsytri2x.f"
		for (j = i__; j <= i__2; ++j) {
#line 360 "dsytri2x.f"
		    a[cut + i__ + (cut + j) * a_dim1] += work[u11 + i__ + j * 
			    work_dim1];
#line 361 "dsytri2x.f"
		}
#line 362 "dsytri2x.f"
	    }

/*        U01 =  U00**T*invD0*U01 */

#line 366 "dsytri2x.f"
	    i__1 = *n + *nb + 1;
#line 366 "dsytri2x.f"
	    dtrmm_("L", uplo, "T", "U", &cut, &nnb, &c_b11, &a[a_offset], lda,
		     &work[work_offset], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);

/*        Update U01 */

#line 372 "dsytri2x.f"
	    i__1 = cut;
#line 372 "dsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 373 "dsytri2x.f"
		i__2 = nnb;
#line 373 "dsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 374 "dsytri2x.f"
		    a[i__ + (cut + j) * a_dim1] = work[i__ + j * work_dim1];
#line 375 "dsytri2x.f"
		}
#line 376 "dsytri2x.f"
	    }

/*      Next Block */

#line 380 "dsytri2x.f"
	}

/*        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */

#line 384 "dsytri2x.f"
	i__ = 1;
#line 385 "dsytri2x.f"
	while(i__ <= *n) {
#line 386 "dsytri2x.f"
	    if (ipiv[i__] > 0) {
#line 387 "dsytri2x.f"
		ip = ipiv[i__];
#line 388 "dsytri2x.f"
		if (i__ < ip) {
#line 388 "dsytri2x.f"
		    dsyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 388 "dsytri2x.f"
		}
#line 389 "dsytri2x.f"
		if (i__ > ip) {
#line 389 "dsytri2x.f"
		    dsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 389 "dsytri2x.f"
		}
#line 390 "dsytri2x.f"
	    } else {
#line 391 "dsytri2x.f"
		ip = -ipiv[i__];
#line 392 "dsytri2x.f"
		++i__;
#line 393 "dsytri2x.f"
		if (i__ - 1 < ip) {
#line 393 "dsytri2x.f"
		    i__1 = i__ - 1;
#line 393 "dsytri2x.f"
		    dsyswapr_(uplo, n, &a[a_offset], lda, &i__1, &ip, (ftnlen)
			    1);
#line 393 "dsytri2x.f"
		}
#line 395 "dsytri2x.f"
		if (i__ - 1 > ip) {
#line 395 "dsytri2x.f"
		    i__1 = i__ - 1;
#line 395 "dsytri2x.f"
		    dsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__1, (ftnlen)
			    1);
#line 395 "dsytri2x.f"
		}
#line 397 "dsytri2x.f"
	    }
#line 398 "dsytri2x.f"
	    ++i__;
#line 399 "dsytri2x.f"
	}
#line 400 "dsytri2x.f"
    } else {

/*        LOWER... */

/*        invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */

#line 406 "dsytri2x.f"
	dtrtri_(uplo, "U", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);

/*       inv(D) and inv(D)*inv(U) */

#line 410 "dsytri2x.f"
	k = *n;
#line 411 "dsytri2x.f"
	while(k >= 1) {
#line 412 "dsytri2x.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal NNB */
#line 414 "dsytri2x.f"
		work[k + invd * work_dim1] = 1. / a[k + k * a_dim1];
#line 415 "dsytri2x.f"
		work[k + (invd + 1) * work_dim1] = 0.;
#line 416 "dsytri2x.f"
		--k;
#line 417 "dsytri2x.f"
	    } else {
/*           2 x 2 diagonal NNB */
#line 419 "dsytri2x.f"
		t = work[k - 1 + work_dim1];
#line 420 "dsytri2x.f"
		ak = a[k - 1 + (k - 1) * a_dim1] / t;
#line 421 "dsytri2x.f"
		akp1 = a[k + k * a_dim1] / t;
#line 422 "dsytri2x.f"
		akkp1 = work[k - 1 + work_dim1] / t;
#line 423 "dsytri2x.f"
		d__ = t * (ak * akp1 - 1.);
#line 424 "dsytri2x.f"
		work[k - 1 + invd * work_dim1] = akp1 / d__;
#line 425 "dsytri2x.f"
		work[k + invd * work_dim1] = ak / d__;
#line 426 "dsytri2x.f"
		work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
#line 427 "dsytri2x.f"
		work[k - 1 + (invd + 1) * work_dim1] = -akkp1 / d__;
#line 428 "dsytri2x.f"
		k += -2;
#line 429 "dsytri2x.f"
	    }
#line 430 "dsytri2x.f"
	}

/*       inv(U**T) = (inv(U))**T */

/*       inv(U**T)*inv(D)*inv(U) */

#line 436 "dsytri2x.f"
	cut = 0;
#line 437 "dsytri2x.f"
	while(cut < *n) {
#line 438 "dsytri2x.f"
	    nnb = *nb;
#line 439 "dsytri2x.f"
	    if (cut + nnb > *n) {
#line 440 "dsytri2x.f"
		nnb = *n - cut;
#line 441 "dsytri2x.f"
	    } else {
#line 442 "dsytri2x.f"
		count = 0;
/*             count negative elements, */
#line 444 "dsytri2x.f"
		i__1 = cut + nnb;
#line 444 "dsytri2x.f"
		for (i__ = cut + 1; i__ <= i__1; ++i__) {
#line 445 "dsytri2x.f"
		    if (ipiv[i__] < 0) {
#line 445 "dsytri2x.f"
			++count;
#line 445 "dsytri2x.f"
		    }
#line 446 "dsytri2x.f"
		}
/*             need a even number for a clear cut */
#line 448 "dsytri2x.f"
		if (count % 2 == 1) {
#line 448 "dsytri2x.f"
		    ++nnb;
#line 448 "dsytri2x.f"
		}
#line 449 "dsytri2x.f"
	    }
/*     L21 Block */
#line 451 "dsytri2x.f"
	    i__1 = *n - cut - nnb;
#line 451 "dsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 452 "dsytri2x.f"
		i__2 = nnb;
#line 452 "dsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 453 "dsytri2x.f"
		    work[i__ + j * work_dim1] = a[cut + nnb + i__ + (cut + j) 
			    * a_dim1];
#line 454 "dsytri2x.f"
		}
#line 455 "dsytri2x.f"
	    }
/*     L11 Block */
#line 457 "dsytri2x.f"
	    i__1 = nnb;
#line 457 "dsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 458 "dsytri2x.f"
		work[u11 + i__ + i__ * work_dim1] = 1.;
#line 459 "dsytri2x.f"
		i__2 = nnb;
#line 459 "dsytri2x.f"
		for (j = i__ + 1; j <= i__2; ++j) {
#line 460 "dsytri2x.f"
		    work[u11 + i__ + j * work_dim1] = 0.;
#line 461 "dsytri2x.f"
		}
#line 462 "dsytri2x.f"
		i__2 = i__ - 1;
#line 462 "dsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 463 "dsytri2x.f"
		    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) 
			    * a_dim1];
#line 464 "dsytri2x.f"
		}
#line 465 "dsytri2x.f"
	    }

/*          invD*L21 */

#line 469 "dsytri2x.f"
	    i__ = *n - cut - nnb;
#line 470 "dsytri2x.f"
	    while(i__ >= 1) {
#line 471 "dsytri2x.f"
		if (ipiv[cut + nnb + i__] > 0) {
#line 472 "dsytri2x.f"
		    i__1 = nnb;
#line 472 "dsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 473 "dsytri2x.f"
			work[i__ + j * work_dim1] = work[cut + nnb + i__ + 
				invd * work_dim1] * work[i__ + j * work_dim1];
#line 474 "dsytri2x.f"
		    }
#line 475 "dsytri2x.f"
		    --i__;
#line 476 "dsytri2x.f"
		} else {
#line 477 "dsytri2x.f"
		    i__1 = nnb;
#line 477 "dsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 478 "dsytri2x.f"
			u01_i_j__ = work[i__ + j * work_dim1];
#line 479 "dsytri2x.f"
			u01_ip1_j__ = work[i__ - 1 + j * work_dim1];
#line 480 "dsytri2x.f"
			work[i__ + j * work_dim1] = work[cut + nnb + i__ + 
				invd * work_dim1] * u01_i_j__ + work[cut + 
				nnb + i__ + (invd + 1) * work_dim1] * 
				u01_ip1_j__;
#line 482 "dsytri2x.f"
			work[i__ - 1 + j * work_dim1] = work[cut + nnb + i__ 
				- 1 + (invd + 1) * work_dim1] * u01_i_j__ + 
				work[cut + nnb + i__ - 1 + invd * work_dim1] *
				 u01_ip1_j__;
#line 484 "dsytri2x.f"
		    }
#line 485 "dsytri2x.f"
		    i__ += -2;
#line 486 "dsytri2x.f"
		}
#line 487 "dsytri2x.f"
	    }

/*        invD1*L11 */

#line 491 "dsytri2x.f"
	    i__ = nnb;
#line 492 "dsytri2x.f"
	    while(i__ >= 1) {
#line 493 "dsytri2x.f"
		if (ipiv[cut + i__] > 0) {
#line 494 "dsytri2x.f"
		    i__1 = nnb;
#line 494 "dsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 495 "dsytri2x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1];
#line 496 "dsytri2x.f"
		    }
#line 497 "dsytri2x.f"
		    --i__;
#line 498 "dsytri2x.f"
		} else {
#line 499 "dsytri2x.f"
		    i__1 = nnb;
#line 499 "dsytri2x.f"
		    for (j = 1; j <= i__1; ++j) {
#line 500 "dsytri2x.f"
			u11_i_j__ = work[u11 + i__ + j * work_dim1];
#line 501 "dsytri2x.f"
			u11_ip1_j__ = work[u11 + i__ - 1 + j * work_dim1];
#line 502 "dsytri2x.f"
			work[u11 + i__ + j * work_dim1] = work[cut + i__ + 
				invd * work_dim1] * work[u11 + i__ + j * 
				work_dim1] + work[cut + i__ + (invd + 1) * 
				work_dim1] * u11_ip1_j__;
#line 504 "dsytri2x.f"
			work[u11 + i__ - 1 + j * work_dim1] = work[cut + i__ 
				- 1 + (invd + 1) * work_dim1] * u11_i_j__ + 
				work[cut + i__ - 1 + invd * work_dim1] * 
				u11_ip1_j__;
#line 506 "dsytri2x.f"
		    }
#line 507 "dsytri2x.f"
		    i__ += -2;
#line 508 "dsytri2x.f"
		}
#line 509 "dsytri2x.f"
	    }

/*       L11**T*invD1*L11->L11 */

#line 513 "dsytri2x.f"
	    i__1 = *n + *nb + 1;
#line 513 "dsytri2x.f"
	    dtrmm_("L", uplo, "T", "U", &nnb, &nnb, &c_b11, &a[cut + 1 + (cut 
		    + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 517 "dsytri2x.f"
	    i__1 = nnb;
#line 517 "dsytri2x.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 518 "dsytri2x.f"
		i__2 = i__;
#line 518 "dsytri2x.f"
		for (j = 1; j <= i__2; ++j) {
#line 519 "dsytri2x.f"
		    a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * 
			    work_dim1];
#line 520 "dsytri2x.f"
		}
#line 521 "dsytri2x.f"
	    }

#line 523 "dsytri2x.f"
	    if (cut + nnb < *n) {

/*          L21**T*invD2*L21->A(CUT+I,CUT+J) */

#line 527 "dsytri2x.f"
		i__1 = *n - nnb - cut;
#line 527 "dsytri2x.f"
		i__2 = *n + *nb + 1;
#line 527 "dsytri2x.f"
		i__3 = *n + *nb + 1;
#line 527 "dsytri2x.f"
		dgemm_("T", "N", &nnb, &nnb, &i__1, &c_b11, &a[cut + nnb + 1 
			+ (cut + 1) * a_dim1], lda, &work[work_offset], &i__2,
			 &c_b15, &work[u11 + 1 + work_dim1], &i__3, (ftnlen)1,
			 (ftnlen)1);

/*        L11 =  L11**T*invD1*L11 + U01**T*invD*U01 */

#line 533 "dsytri2x.f"
		i__1 = nnb;
#line 533 "dsytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 534 "dsytri2x.f"
		    i__2 = i__;
#line 534 "dsytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 535 "dsytri2x.f"
			a[cut + i__ + (cut + j) * a_dim1] += work[u11 + i__ + 
				j * work_dim1];
#line 536 "dsytri2x.f"
		    }
#line 537 "dsytri2x.f"
		}

/*        L01 =  L22**T*invD2*L21 */

#line 541 "dsytri2x.f"
		i__1 = *n - nnb - cut;
#line 541 "dsytri2x.f"
		i__2 = *n + *nb + 1;
#line 541 "dsytri2x.f"
		dtrmm_("L", uplo, "T", "U", &i__1, &nnb, &c_b11, &a[cut + nnb 
			+ 1 + (cut + nnb + 1) * a_dim1], lda, &work[
			work_offset], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);

/*      Update L21 */

#line 546 "dsytri2x.f"
		i__1 = *n - cut - nnb;
#line 546 "dsytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 547 "dsytri2x.f"
		    i__2 = nnb;
#line 547 "dsytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 548 "dsytri2x.f"
			a[cut + nnb + i__ + (cut + j) * a_dim1] = work[i__ + 
				j * work_dim1];
#line 549 "dsytri2x.f"
		    }
#line 550 "dsytri2x.f"
		}
#line 552 "dsytri2x.f"
	    } else {

/*        L11 =  L11**T*invD1*L11 */

#line 556 "dsytri2x.f"
		i__1 = nnb;
#line 556 "dsytri2x.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 557 "dsytri2x.f"
		    i__2 = i__;
#line 557 "dsytri2x.f"
		    for (j = 1; j <= i__2; ++j) {
#line 558 "dsytri2x.f"
			a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + 
				j * work_dim1];
#line 559 "dsytri2x.f"
		    }
#line 560 "dsytri2x.f"
		}
#line 561 "dsytri2x.f"
	    }

/*      Next Block */

#line 565 "dsytri2x.f"
	    cut += nnb;
#line 566 "dsytri2x.f"
	}

/*        Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */

#line 570 "dsytri2x.f"
	i__ = *n;
#line 571 "dsytri2x.f"
	while(i__ >= 1) {
#line 572 "dsytri2x.f"
	    if (ipiv[i__] > 0) {
#line 573 "dsytri2x.f"
		ip = ipiv[i__];
#line 574 "dsytri2x.f"
		if (i__ < ip) {
#line 574 "dsytri2x.f"
		    dsyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 574 "dsytri2x.f"
		}
#line 575 "dsytri2x.f"
		if (i__ > ip) {
#line 575 "dsytri2x.f"
		    dsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 575 "dsytri2x.f"
		}
#line 576 "dsytri2x.f"
	    } else {
#line 577 "dsytri2x.f"
		ip = -ipiv[i__];
#line 578 "dsytri2x.f"
		if (i__ < ip) {
#line 578 "dsytri2x.f"
		    dsyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip, (ftnlen)
			    1);
#line 578 "dsytri2x.f"
		}
#line 579 "dsytri2x.f"
		if (i__ > ip) {
#line 579 "dsytri2x.f"
		    dsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__, (ftnlen)
			    1);
#line 579 "dsytri2x.f"
		}
#line 580 "dsytri2x.f"
		--i__;
#line 581 "dsytri2x.f"
	    }
#line 582 "dsytri2x.f"
	    --i__;
#line 583 "dsytri2x.f"
	}
#line 584 "dsytri2x.f"
    }

#line 586 "dsytri2x.f"
    return 0;

/*     End of DSYTRI2X */

} /* dsytri2x_ */

