#line 1 "zsytri_rook.f"
/* zsytri_rook.f -- translated by f2c (version 20100827).
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

#line 1 "zsytri_rook.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b ZSYTRI_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRI_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYTRI_ROOK computes the inverse of a complex symmetric */
/* > matrix A using the factorization A = U*D*U**T or A = L*D*L**T */
/* > computed by ZSYTRF_ROOK. */
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
/* >          On entry, the block diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by ZSYTRF_ROOK. */
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
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZSYTRF_ROOK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) */
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

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >   December 2016, Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zsytri_rook__(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex d__;
    static integer k;
    static doublecomplex t, ak;
    static integer kp;
    static doublecomplex akp1, temp, akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID zdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zsymv_(char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

#line 173 "zsytri_rook.f"
    /* Parameter adjustments */
#line 173 "zsytri_rook.f"
    a_dim1 = *lda;
#line 173 "zsytri_rook.f"
    a_offset = 1 + a_dim1;
#line 173 "zsytri_rook.f"
    a -= a_offset;
#line 173 "zsytri_rook.f"
    --ipiv;
#line 173 "zsytri_rook.f"
    --work;
#line 173 "zsytri_rook.f"

#line 173 "zsytri_rook.f"
    /* Function Body */
#line 173 "zsytri_rook.f"
    *info = 0;
#line 174 "zsytri_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 175 "zsytri_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 176 "zsytri_rook.f"
	*info = -1;
#line 177 "zsytri_rook.f"
    } else if (*n < 0) {
#line 178 "zsytri_rook.f"
	*info = -2;
#line 179 "zsytri_rook.f"
    } else if (*lda < max(1,*n)) {
#line 180 "zsytri_rook.f"
	*info = -4;
#line 181 "zsytri_rook.f"
    }
#line 182 "zsytri_rook.f"
    if (*info != 0) {
#line 183 "zsytri_rook.f"
	i__1 = -(*info);
#line 183 "zsytri_rook.f"
	xerbla_("ZSYTRI_ROOK", &i__1, (ftnlen)11);
#line 184 "zsytri_rook.f"
	return 0;
#line 185 "zsytri_rook.f"
    }

/*     Quick return if possible */

#line 189 "zsytri_rook.f"
    if (*n == 0) {
#line 189 "zsytri_rook.f"
	return 0;
#line 189 "zsytri_rook.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 194 "zsytri_rook.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 198 "zsytri_rook.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 199 "zsytri_rook.f"
	    i__1 = *info + *info * a_dim1;
#line 199 "zsytri_rook.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 199 "zsytri_rook.f"
		return 0;
#line 199 "zsytri_rook.f"
	    }
#line 201 "zsytri_rook.f"
/* L10: */
#line 201 "zsytri_rook.f"
	}
#line 202 "zsytri_rook.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 206 "zsytri_rook.f"
	i__1 = *n;
#line 206 "zsytri_rook.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 207 "zsytri_rook.f"
	    i__2 = *info + *info * a_dim1;
#line 207 "zsytri_rook.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 207 "zsytri_rook.f"
		return 0;
#line 207 "zsytri_rook.f"
	    }
#line 209 "zsytri_rook.f"
/* L20: */
#line 209 "zsytri_rook.f"
	}
#line 210 "zsytri_rook.f"
    }
#line 211 "zsytri_rook.f"
    *info = 0;

#line 213 "zsytri_rook.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 220 "zsytri_rook.f"
	k = 1;
#line 221 "zsytri_rook.f"
L30:

/*        If K > N, exit from loop. */

#line 225 "zsytri_rook.f"
	if (k > *n) {
#line 225 "zsytri_rook.f"
	    goto L40;
#line 225 "zsytri_rook.f"
	}

#line 228 "zsytri_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 234 "zsytri_rook.f"
	    i__1 = k + k * a_dim1;
#line 234 "zsytri_rook.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 234 "zsytri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute column K of the inverse. */

#line 238 "zsytri_rook.f"
	    if (k > 1) {
#line 239 "zsytri_rook.f"
		i__1 = k - 1;
#line 239 "zsytri_rook.f"
		zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 240 "zsytri_rook.f"
		i__1 = k - 1;
#line 240 "zsytri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 240 "zsytri_rook.f"
		zsymv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 242 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 242 "zsytri_rook.f"
		i__2 = k + k * a_dim1;
#line 242 "zsytri_rook.f"
		i__3 = k - 1;
#line 242 "zsytri_rook.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 242 "zsytri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 242 "zsytri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 244 "zsytri_rook.f"
	    }
#line 245 "zsytri_rook.f"
	    kstep = 1;
#line 246 "zsytri_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 252 "zsytri_rook.f"
	    i__1 = k + (k + 1) * a_dim1;
#line 252 "zsytri_rook.f"
	    t.r = a[i__1].r, t.i = a[i__1].i;
#line 253 "zsytri_rook.f"
	    z_div(&z__1, &a[k + k * a_dim1], &t);
#line 253 "zsytri_rook.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 254 "zsytri_rook.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &t);
#line 254 "zsytri_rook.f"
	    akp1.r = z__1.r, akp1.i = z__1.i;
#line 255 "zsytri_rook.f"
	    z_div(&z__1, &a[k + (k + 1) * a_dim1], &t);
#line 255 "zsytri_rook.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 256 "zsytri_rook.f"
	    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + 
		    ak.i * akp1.r;
#line 256 "zsytri_rook.f"
	    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 256 "zsytri_rook.f"
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
#line 256 "zsytri_rook.f"
	    d__.r = z__1.r, d__.i = z__1.i;
#line 257 "zsytri_rook.f"
	    i__1 = k + k * a_dim1;
#line 257 "zsytri_rook.f"
	    z_div(&z__1, &akp1, &d__);
#line 257 "zsytri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 258 "zsytri_rook.f"
	    i__1 = k + 1 + (k + 1) * a_dim1;
#line 258 "zsytri_rook.f"
	    z_div(&z__1, &ak, &d__);
#line 258 "zsytri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 259 "zsytri_rook.f"
	    i__1 = k + (k + 1) * a_dim1;
#line 259 "zsytri_rook.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 259 "zsytri_rook.f"
	    z_div(&z__1, &z__2, &d__);
#line 259 "zsytri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute columns K and K+1 of the inverse. */

#line 263 "zsytri_rook.f"
	    if (k > 1) {
#line 264 "zsytri_rook.f"
		i__1 = k - 1;
#line 264 "zsytri_rook.f"
		zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 265 "zsytri_rook.f"
		i__1 = k - 1;
#line 265 "zsytri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 265 "zsytri_rook.f"
		zsymv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 267 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 267 "zsytri_rook.f"
		i__2 = k + k * a_dim1;
#line 267 "zsytri_rook.f"
		i__3 = k - 1;
#line 267 "zsytri_rook.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 267 "zsytri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 267 "zsytri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 269 "zsytri_rook.f"
		i__1 = k + (k + 1) * a_dim1;
#line 269 "zsytri_rook.f"
		i__2 = k + (k + 1) * a_dim1;
#line 269 "zsytri_rook.f"
		i__3 = k - 1;
#line 269 "zsytri_rook.f"
		zdotu_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * 
			a_dim1 + 1], &c__1);
#line 269 "zsytri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 269 "zsytri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 271 "zsytri_rook.f"
		i__1 = k - 1;
#line 271 "zsytri_rook.f"
		zcopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &
			c__1);
#line 272 "zsytri_rook.f"
		i__1 = k - 1;
#line 272 "zsytri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 272 "zsytri_rook.f"
		zsymv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[(k + 1) * a_dim1 + 1], &c__1, (ftnlen)1);
#line 274 "zsytri_rook.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 274 "zsytri_rook.f"
		i__2 = k + 1 + (k + 1) * a_dim1;
#line 274 "zsytri_rook.f"
		i__3 = k - 1;
#line 274 "zsytri_rook.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1]
			, &c__1);
#line 274 "zsytri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 274 "zsytri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 276 "zsytri_rook.f"
	    }
#line 277 "zsytri_rook.f"
	    kstep = 2;
#line 278 "zsytri_rook.f"
	}

#line 280 "zsytri_rook.f"
	if (kstep == 1) {

/*           Interchange rows and columns K and IPIV(K) in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 285 "zsytri_rook.f"
	    kp = ipiv[k];
#line 286 "zsytri_rook.f"
	    if (kp != k) {
#line 287 "zsytri_rook.f"
		if (kp > 1) {
#line 287 "zsytri_rook.f"
		    i__1 = kp - 1;
#line 287 "zsytri_rook.f"
		    zswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 287 "zsytri_rook.f"
		}
#line 289 "zsytri_rook.f"
		i__1 = k - kp - 1;
#line 289 "zsytri_rook.f"
		zswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);
#line 290 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 290 "zsytri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 291 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 291 "zsytri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 291 "zsytri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 292 "zsytri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 292 "zsytri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 293 "zsytri_rook.f"
	    }
#line 294 "zsytri_rook.f"
	} else {

/*           Interchange rows and columns K and K+1 with -IPIV(K) and */
/*           -IPIV(K+1)in the leading submatrix A(1:k+1,1:k+1) */

#line 299 "zsytri_rook.f"
	    kp = -ipiv[k];
#line 300 "zsytri_rook.f"
	    if (kp != k) {
#line 301 "zsytri_rook.f"
		if (kp > 1) {
#line 301 "zsytri_rook.f"
		    i__1 = kp - 1;
#line 301 "zsytri_rook.f"
		    zswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 301 "zsytri_rook.f"
		}
#line 303 "zsytri_rook.f"
		i__1 = k - kp - 1;
#line 303 "zsytri_rook.f"
		zswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);

#line 305 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 305 "zsytri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 306 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 306 "zsytri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 306 "zsytri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 307 "zsytri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 307 "zsytri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 308 "zsytri_rook.f"
		i__1 = k + (k + 1) * a_dim1;
#line 308 "zsytri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 309 "zsytri_rook.f"
		i__1 = k + (k + 1) * a_dim1;
#line 309 "zsytri_rook.f"
		i__2 = kp + (k + 1) * a_dim1;
#line 309 "zsytri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 310 "zsytri_rook.f"
		i__1 = kp + (k + 1) * a_dim1;
#line 310 "zsytri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 311 "zsytri_rook.f"
	    }

#line 313 "zsytri_rook.f"
	    ++k;
#line 314 "zsytri_rook.f"
	    kp = -ipiv[k];
#line 315 "zsytri_rook.f"
	    if (kp != k) {
#line 316 "zsytri_rook.f"
		if (kp > 1) {
#line 316 "zsytri_rook.f"
		    i__1 = kp - 1;
#line 316 "zsytri_rook.f"
		    zswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 316 "zsytri_rook.f"
		}
#line 318 "zsytri_rook.f"
		i__1 = k - kp - 1;
#line 318 "zsytri_rook.f"
		zswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);
#line 319 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 319 "zsytri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 320 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 320 "zsytri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 320 "zsytri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 321 "zsytri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 321 "zsytri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 322 "zsytri_rook.f"
	    }
#line 323 "zsytri_rook.f"
	}

#line 325 "zsytri_rook.f"
	++k;
#line 326 "zsytri_rook.f"
	goto L30;
#line 327 "zsytri_rook.f"
L40:

#line 329 "zsytri_rook.f"
	;
#line 329 "zsytri_rook.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 336 "zsytri_rook.f"
	k = *n;
#line 337 "zsytri_rook.f"
L50:

/*        If K < 1, exit from loop. */

#line 341 "zsytri_rook.f"
	if (k < 1) {
#line 341 "zsytri_rook.f"
	    goto L60;
#line 341 "zsytri_rook.f"
	}

#line 344 "zsytri_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 350 "zsytri_rook.f"
	    i__1 = k + k * a_dim1;
#line 350 "zsytri_rook.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 350 "zsytri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute column K of the inverse. */

#line 354 "zsytri_rook.f"
	    if (k < *n) {
#line 355 "zsytri_rook.f"
		i__1 = *n - k;
#line 355 "zsytri_rook.f"
		zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 356 "zsytri_rook.f"
		i__1 = *n - k;
#line 356 "zsytri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 356 "zsytri_rook.f"
		zsymv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1,
			 (ftnlen)1);
#line 358 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 358 "zsytri_rook.f"
		i__2 = k + k * a_dim1;
#line 358 "zsytri_rook.f"
		i__3 = *n - k;
#line 358 "zsytri_rook.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], 
			&c__1);
#line 358 "zsytri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 358 "zsytri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 360 "zsytri_rook.f"
	    }
#line 361 "zsytri_rook.f"
	    kstep = 1;
#line 362 "zsytri_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 368 "zsytri_rook.f"
	    i__1 = k + (k - 1) * a_dim1;
#line 368 "zsytri_rook.f"
	    t.r = a[i__1].r, t.i = a[i__1].i;
#line 369 "zsytri_rook.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &t);
#line 369 "zsytri_rook.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 370 "zsytri_rook.f"
	    z_div(&z__1, &a[k + k * a_dim1], &t);
#line 370 "zsytri_rook.f"
	    akp1.r = z__1.r, akp1.i = z__1.i;
#line 371 "zsytri_rook.f"
	    z_div(&z__1, &a[k + (k - 1) * a_dim1], &t);
#line 371 "zsytri_rook.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 372 "zsytri_rook.f"
	    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + 
		    ak.i * akp1.r;
#line 372 "zsytri_rook.f"
	    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 372 "zsytri_rook.f"
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
#line 372 "zsytri_rook.f"
	    d__.r = z__1.r, d__.i = z__1.i;
#line 373 "zsytri_rook.f"
	    i__1 = k - 1 + (k - 1) * a_dim1;
#line 373 "zsytri_rook.f"
	    z_div(&z__1, &akp1, &d__);
#line 373 "zsytri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 374 "zsytri_rook.f"
	    i__1 = k + k * a_dim1;
#line 374 "zsytri_rook.f"
	    z_div(&z__1, &ak, &d__);
#line 374 "zsytri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 375 "zsytri_rook.f"
	    i__1 = k + (k - 1) * a_dim1;
#line 375 "zsytri_rook.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 375 "zsytri_rook.f"
	    z_div(&z__1, &z__2, &d__);
#line 375 "zsytri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute columns K-1 and K of the inverse. */

#line 379 "zsytri_rook.f"
	    if (k < *n) {
#line 380 "zsytri_rook.f"
		i__1 = *n - k;
#line 380 "zsytri_rook.f"
		zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 381 "zsytri_rook.f"
		i__1 = *n - k;
#line 381 "zsytri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 381 "zsytri_rook.f"
		zsymv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1,
			 (ftnlen)1);
#line 383 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 383 "zsytri_rook.f"
		i__2 = k + k * a_dim1;
#line 383 "zsytri_rook.f"
		i__3 = *n - k;
#line 383 "zsytri_rook.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], 
			&c__1);
#line 383 "zsytri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 383 "zsytri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 385 "zsytri_rook.f"
		i__1 = k + (k - 1) * a_dim1;
#line 385 "zsytri_rook.f"
		i__2 = k + (k - 1) * a_dim1;
#line 385 "zsytri_rook.f"
		i__3 = *n - k;
#line 385 "zsytri_rook.f"
		zdotu_(&z__2, &i__3, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 
			+ (k - 1) * a_dim1], &c__1);
#line 385 "zsytri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 385 "zsytri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 388 "zsytri_rook.f"
		i__1 = *n - k;
#line 388 "zsytri_rook.f"
		zcopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &
			c__1);
#line 389 "zsytri_rook.f"
		i__1 = *n - k;
#line 389 "zsytri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 389 "zsytri_rook.f"
		zsymv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + (k - 1) * a_dim1], 
			&c__1, (ftnlen)1);
#line 391 "zsytri_rook.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 391 "zsytri_rook.f"
		i__2 = k - 1 + (k - 1) * a_dim1;
#line 391 "zsytri_rook.f"
		i__3 = *n - k;
#line 391 "zsytri_rook.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + (k - 1) * 
			a_dim1], &c__1);
#line 391 "zsytri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 391 "zsytri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 393 "zsytri_rook.f"
	    }
#line 394 "zsytri_rook.f"
	    kstep = 2;
#line 395 "zsytri_rook.f"
	}

#line 397 "zsytri_rook.f"
	if (kstep == 1) {

/*           Interchange rows and columns K and IPIV(K) in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 402 "zsytri_rook.f"
	    kp = ipiv[k];
#line 403 "zsytri_rook.f"
	    if (kp != k) {
#line 404 "zsytri_rook.f"
		if (kp < *n) {
#line 404 "zsytri_rook.f"
		    i__1 = *n - kp;
#line 404 "zsytri_rook.f"
		    zswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 404 "zsytri_rook.f"
		}
#line 406 "zsytri_rook.f"
		i__1 = kp - k - 1;
#line 406 "zsytri_rook.f"
		zswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);
#line 407 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 407 "zsytri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 408 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 408 "zsytri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 408 "zsytri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 409 "zsytri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 409 "zsytri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 410 "zsytri_rook.f"
	    }
#line 411 "zsytri_rook.f"
	} else {

/*           Interchange rows and columns K and K-1 with -IPIV(K) and */
/*           -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n) */

#line 416 "zsytri_rook.f"
	    kp = -ipiv[k];
#line 417 "zsytri_rook.f"
	    if (kp != k) {
#line 418 "zsytri_rook.f"
		if (kp < *n) {
#line 418 "zsytri_rook.f"
		    i__1 = *n - kp;
#line 418 "zsytri_rook.f"
		    zswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 418 "zsytri_rook.f"
		}
#line 420 "zsytri_rook.f"
		i__1 = kp - k - 1;
#line 420 "zsytri_rook.f"
		zswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);

#line 422 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 422 "zsytri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 423 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 423 "zsytri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 423 "zsytri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 424 "zsytri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 424 "zsytri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 425 "zsytri_rook.f"
		i__1 = k + (k - 1) * a_dim1;
#line 425 "zsytri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 426 "zsytri_rook.f"
		i__1 = k + (k - 1) * a_dim1;
#line 426 "zsytri_rook.f"
		i__2 = kp + (k - 1) * a_dim1;
#line 426 "zsytri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 427 "zsytri_rook.f"
		i__1 = kp + (k - 1) * a_dim1;
#line 427 "zsytri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 428 "zsytri_rook.f"
	    }

#line 430 "zsytri_rook.f"
	    --k;
#line 431 "zsytri_rook.f"
	    kp = -ipiv[k];
#line 432 "zsytri_rook.f"
	    if (kp != k) {
#line 433 "zsytri_rook.f"
		if (kp < *n) {
#line 433 "zsytri_rook.f"
		    i__1 = *n - kp;
#line 433 "zsytri_rook.f"
		    zswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 433 "zsytri_rook.f"
		}
#line 435 "zsytri_rook.f"
		i__1 = kp - k - 1;
#line 435 "zsytri_rook.f"
		zswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);
#line 436 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 436 "zsytri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 437 "zsytri_rook.f"
		i__1 = k + k * a_dim1;
#line 437 "zsytri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 437 "zsytri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 438 "zsytri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 438 "zsytri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 439 "zsytri_rook.f"
	    }
#line 440 "zsytri_rook.f"
	}

#line 442 "zsytri_rook.f"
	--k;
#line 443 "zsytri_rook.f"
	goto L50;
#line 444 "zsytri_rook.f"
L60:
#line 445 "zsytri_rook.f"
	;
#line 445 "zsytri_rook.f"
    }

#line 447 "zsytri_rook.f"
    return 0;

/*     End of ZSYTRI_ROOK */

} /* zsytri_rook__ */

