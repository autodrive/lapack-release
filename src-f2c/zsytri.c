#line 1 "zsytri.f"
/* zsytri.f -- translated by f2c (version 20100827).
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

#line 1 "zsytri.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b ZSYTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO ) */

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
/* > ZSYTRI computes the inverse of a complex symmetric indefinite matrix */
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
/* >          On entry, the block diagonal matrix D and the multipliers */
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
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
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

/* > \date November 2011 */

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsytri_(char *uplo, integer *n, doublecomplex *a, 
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


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 158 "zsytri.f"
    /* Parameter adjustments */
#line 158 "zsytri.f"
    a_dim1 = *lda;
#line 158 "zsytri.f"
    a_offset = 1 + a_dim1;
#line 158 "zsytri.f"
    a -= a_offset;
#line 158 "zsytri.f"
    --ipiv;
#line 158 "zsytri.f"
    --work;
#line 158 "zsytri.f"

#line 158 "zsytri.f"
    /* Function Body */
#line 158 "zsytri.f"
    *info = 0;
#line 159 "zsytri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 160 "zsytri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 161 "zsytri.f"
	*info = -1;
#line 162 "zsytri.f"
    } else if (*n < 0) {
#line 163 "zsytri.f"
	*info = -2;
#line 164 "zsytri.f"
    } else if (*lda < max(1,*n)) {
#line 165 "zsytri.f"
	*info = -4;
#line 166 "zsytri.f"
    }
#line 167 "zsytri.f"
    if (*info != 0) {
#line 168 "zsytri.f"
	i__1 = -(*info);
#line 168 "zsytri.f"
	xerbla_("ZSYTRI", &i__1, (ftnlen)6);
#line 169 "zsytri.f"
	return 0;
#line 170 "zsytri.f"
    }

/*     Quick return if possible */

#line 174 "zsytri.f"
    if (*n == 0) {
#line 174 "zsytri.f"
	return 0;
#line 174 "zsytri.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 179 "zsytri.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 183 "zsytri.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 184 "zsytri.f"
	    i__1 = *info + *info * a_dim1;
#line 184 "zsytri.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 184 "zsytri.f"
		return 0;
#line 184 "zsytri.f"
	    }
#line 186 "zsytri.f"
/* L10: */
#line 186 "zsytri.f"
	}
#line 187 "zsytri.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 191 "zsytri.f"
	i__1 = *n;
#line 191 "zsytri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 192 "zsytri.f"
	    i__2 = *info + *info * a_dim1;
#line 192 "zsytri.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 192 "zsytri.f"
		return 0;
#line 192 "zsytri.f"
	    }
#line 194 "zsytri.f"
/* L20: */
#line 194 "zsytri.f"
	}
#line 195 "zsytri.f"
    }
#line 196 "zsytri.f"
    *info = 0;

#line 198 "zsytri.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 205 "zsytri.f"
	k = 1;
#line 206 "zsytri.f"
L30:

/*        If K > N, exit from loop. */

#line 210 "zsytri.f"
	if (k > *n) {
#line 210 "zsytri.f"
	    goto L40;
#line 210 "zsytri.f"
	}

#line 213 "zsytri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 219 "zsytri.f"
	    i__1 = k + k * a_dim1;
#line 219 "zsytri.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 219 "zsytri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute column K of the inverse. */

#line 223 "zsytri.f"
	    if (k > 1) {
#line 224 "zsytri.f"
		i__1 = k - 1;
#line 224 "zsytri.f"
		zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 225 "zsytri.f"
		i__1 = k - 1;
#line 225 "zsytri.f"
		z__1.r = -1., z__1.i = -0.;
#line 225 "zsytri.f"
		zsymv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 227 "zsytri.f"
		i__1 = k + k * a_dim1;
#line 227 "zsytri.f"
		i__2 = k + k * a_dim1;
#line 227 "zsytri.f"
		i__3 = k - 1;
#line 227 "zsytri.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 227 "zsytri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 227 "zsytri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 229 "zsytri.f"
	    }
#line 230 "zsytri.f"
	    kstep = 1;
#line 231 "zsytri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 237 "zsytri.f"
	    i__1 = k + (k + 1) * a_dim1;
#line 237 "zsytri.f"
	    t.r = a[i__1].r, t.i = a[i__1].i;
#line 238 "zsytri.f"
	    z_div(&z__1, &a[k + k * a_dim1], &t);
#line 238 "zsytri.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 239 "zsytri.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &t);
#line 239 "zsytri.f"
	    akp1.r = z__1.r, akp1.i = z__1.i;
#line 240 "zsytri.f"
	    z_div(&z__1, &a[k + (k + 1) * a_dim1], &t);
#line 240 "zsytri.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 241 "zsytri.f"
	    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + 
		    ak.i * akp1.r;
#line 241 "zsytri.f"
	    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 241 "zsytri.f"
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
#line 241 "zsytri.f"
	    d__.r = z__1.r, d__.i = z__1.i;
#line 242 "zsytri.f"
	    i__1 = k + k * a_dim1;
#line 242 "zsytri.f"
	    z_div(&z__1, &akp1, &d__);
#line 242 "zsytri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 243 "zsytri.f"
	    i__1 = k + 1 + (k + 1) * a_dim1;
#line 243 "zsytri.f"
	    z_div(&z__1, &ak, &d__);
#line 243 "zsytri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 244 "zsytri.f"
	    i__1 = k + (k + 1) * a_dim1;
#line 244 "zsytri.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 244 "zsytri.f"
	    z_div(&z__1, &z__2, &d__);
#line 244 "zsytri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute columns K and K+1 of the inverse. */

#line 248 "zsytri.f"
	    if (k > 1) {
#line 249 "zsytri.f"
		i__1 = k - 1;
#line 249 "zsytri.f"
		zcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 250 "zsytri.f"
		i__1 = k - 1;
#line 250 "zsytri.f"
		z__1.r = -1., z__1.i = -0.;
#line 250 "zsytri.f"
		zsymv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 252 "zsytri.f"
		i__1 = k + k * a_dim1;
#line 252 "zsytri.f"
		i__2 = k + k * a_dim1;
#line 252 "zsytri.f"
		i__3 = k - 1;
#line 252 "zsytri.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 252 "zsytri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 252 "zsytri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 254 "zsytri.f"
		i__1 = k + (k + 1) * a_dim1;
#line 254 "zsytri.f"
		i__2 = k + (k + 1) * a_dim1;
#line 254 "zsytri.f"
		i__3 = k - 1;
#line 254 "zsytri.f"
		zdotu_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * 
			a_dim1 + 1], &c__1);
#line 254 "zsytri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 254 "zsytri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 256 "zsytri.f"
		i__1 = k - 1;
#line 256 "zsytri.f"
		zcopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &
			c__1);
#line 257 "zsytri.f"
		i__1 = k - 1;
#line 257 "zsytri.f"
		z__1.r = -1., z__1.i = -0.;
#line 257 "zsytri.f"
		zsymv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[(k + 1) * a_dim1 + 1], &c__1, (ftnlen)1);
#line 259 "zsytri.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 259 "zsytri.f"
		i__2 = k + 1 + (k + 1) * a_dim1;
#line 259 "zsytri.f"
		i__3 = k - 1;
#line 259 "zsytri.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1]
			, &c__1);
#line 259 "zsytri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 259 "zsytri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 261 "zsytri.f"
	    }
#line 262 "zsytri.f"
	    kstep = 2;
#line 263 "zsytri.f"
	}

#line 265 "zsytri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 266 "zsytri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 271 "zsytri.f"
	    i__1 = kp - 1;
#line 271 "zsytri.f"
	    zswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
		    c__1);
#line 272 "zsytri.f"
	    i__1 = k - kp - 1;
#line 272 "zsytri.f"
	    zswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * 
		    a_dim1], lda);
#line 273 "zsytri.f"
	    i__1 = k + k * a_dim1;
#line 273 "zsytri.f"
	    temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 274 "zsytri.f"
	    i__1 = k + k * a_dim1;
#line 274 "zsytri.f"
	    i__2 = kp + kp * a_dim1;
#line 274 "zsytri.f"
	    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 275 "zsytri.f"
	    i__1 = kp + kp * a_dim1;
#line 275 "zsytri.f"
	    a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 276 "zsytri.f"
	    if (kstep == 2) {
#line 277 "zsytri.f"
		i__1 = k + (k + 1) * a_dim1;
#line 277 "zsytri.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 278 "zsytri.f"
		i__1 = k + (k + 1) * a_dim1;
#line 278 "zsytri.f"
		i__2 = kp + (k + 1) * a_dim1;
#line 278 "zsytri.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 279 "zsytri.f"
		i__1 = kp + (k + 1) * a_dim1;
#line 279 "zsytri.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 280 "zsytri.f"
	    }
#line 281 "zsytri.f"
	}

#line 283 "zsytri.f"
	k += kstep;
#line 284 "zsytri.f"
	goto L30;
#line 285 "zsytri.f"
L40:

#line 287 "zsytri.f"
	;
#line 287 "zsytri.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 294 "zsytri.f"
	k = *n;
#line 295 "zsytri.f"
L50:

/*        If K < 1, exit from loop. */

#line 299 "zsytri.f"
	if (k < 1) {
#line 299 "zsytri.f"
	    goto L60;
#line 299 "zsytri.f"
	}

#line 302 "zsytri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 308 "zsytri.f"
	    i__1 = k + k * a_dim1;
#line 308 "zsytri.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 308 "zsytri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute column K of the inverse. */

#line 312 "zsytri.f"
	    if (k < *n) {
#line 313 "zsytri.f"
		i__1 = *n - k;
#line 313 "zsytri.f"
		zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 314 "zsytri.f"
		i__1 = *n - k;
#line 314 "zsytri.f"
		z__1.r = -1., z__1.i = -0.;
#line 314 "zsytri.f"
		zsymv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1,
			 (ftnlen)1);
#line 316 "zsytri.f"
		i__1 = k + k * a_dim1;
#line 316 "zsytri.f"
		i__2 = k + k * a_dim1;
#line 316 "zsytri.f"
		i__3 = *n - k;
#line 316 "zsytri.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], 
			&c__1);
#line 316 "zsytri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 316 "zsytri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 318 "zsytri.f"
	    }
#line 319 "zsytri.f"
	    kstep = 1;
#line 320 "zsytri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 326 "zsytri.f"
	    i__1 = k + (k - 1) * a_dim1;
#line 326 "zsytri.f"
	    t.r = a[i__1].r, t.i = a[i__1].i;
#line 327 "zsytri.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &t);
#line 327 "zsytri.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 328 "zsytri.f"
	    z_div(&z__1, &a[k + k * a_dim1], &t);
#line 328 "zsytri.f"
	    akp1.r = z__1.r, akp1.i = z__1.i;
#line 329 "zsytri.f"
	    z_div(&z__1, &a[k + (k - 1) * a_dim1], &t);
#line 329 "zsytri.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 330 "zsytri.f"
	    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + 
		    ak.i * akp1.r;
#line 330 "zsytri.f"
	    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 330 "zsytri.f"
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
#line 330 "zsytri.f"
	    d__.r = z__1.r, d__.i = z__1.i;
#line 331 "zsytri.f"
	    i__1 = k - 1 + (k - 1) * a_dim1;
#line 331 "zsytri.f"
	    z_div(&z__1, &akp1, &d__);
#line 331 "zsytri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 332 "zsytri.f"
	    i__1 = k + k * a_dim1;
#line 332 "zsytri.f"
	    z_div(&z__1, &ak, &d__);
#line 332 "zsytri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 333 "zsytri.f"
	    i__1 = k + (k - 1) * a_dim1;
#line 333 "zsytri.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 333 "zsytri.f"
	    z_div(&z__1, &z__2, &d__);
#line 333 "zsytri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute columns K-1 and K of the inverse. */

#line 337 "zsytri.f"
	    if (k < *n) {
#line 338 "zsytri.f"
		i__1 = *n - k;
#line 338 "zsytri.f"
		zcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 339 "zsytri.f"
		i__1 = *n - k;
#line 339 "zsytri.f"
		z__1.r = -1., z__1.i = -0.;
#line 339 "zsytri.f"
		zsymv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1,
			 (ftnlen)1);
#line 341 "zsytri.f"
		i__1 = k + k * a_dim1;
#line 341 "zsytri.f"
		i__2 = k + k * a_dim1;
#line 341 "zsytri.f"
		i__3 = *n - k;
#line 341 "zsytri.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], 
			&c__1);
#line 341 "zsytri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 341 "zsytri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 343 "zsytri.f"
		i__1 = k + (k - 1) * a_dim1;
#line 343 "zsytri.f"
		i__2 = k + (k - 1) * a_dim1;
#line 343 "zsytri.f"
		i__3 = *n - k;
#line 343 "zsytri.f"
		zdotu_(&z__2, &i__3, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 
			+ (k - 1) * a_dim1], &c__1);
#line 343 "zsytri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 343 "zsytri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 346 "zsytri.f"
		i__1 = *n - k;
#line 346 "zsytri.f"
		zcopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &
			c__1);
#line 347 "zsytri.f"
		i__1 = *n - k;
#line 347 "zsytri.f"
		z__1.r = -1., z__1.i = -0.;
#line 347 "zsytri.f"
		zsymv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + (k - 1) * a_dim1], 
			&c__1, (ftnlen)1);
#line 349 "zsytri.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 349 "zsytri.f"
		i__2 = k - 1 + (k - 1) * a_dim1;
#line 349 "zsytri.f"
		i__3 = *n - k;
#line 349 "zsytri.f"
		zdotu_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + (k - 1) * 
			a_dim1], &c__1);
#line 349 "zsytri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 349 "zsytri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 351 "zsytri.f"
	    }
#line 352 "zsytri.f"
	    kstep = 2;
#line 353 "zsytri.f"
	}

#line 355 "zsytri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 356 "zsytri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 361 "zsytri.f"
	    if (kp < *n) {
#line 361 "zsytri.f"
		i__1 = *n - kp;
#line 361 "zsytri.f"
		zswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp *
			 a_dim1], &c__1);
#line 361 "zsytri.f"
	    }
#line 363 "zsytri.f"
	    i__1 = kp - k - 1;
#line 363 "zsytri.f"
	    zswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * 
		    a_dim1], lda);
#line 364 "zsytri.f"
	    i__1 = k + k * a_dim1;
#line 364 "zsytri.f"
	    temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 365 "zsytri.f"
	    i__1 = k + k * a_dim1;
#line 365 "zsytri.f"
	    i__2 = kp + kp * a_dim1;
#line 365 "zsytri.f"
	    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 366 "zsytri.f"
	    i__1 = kp + kp * a_dim1;
#line 366 "zsytri.f"
	    a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 367 "zsytri.f"
	    if (kstep == 2) {
#line 368 "zsytri.f"
		i__1 = k + (k - 1) * a_dim1;
#line 368 "zsytri.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 369 "zsytri.f"
		i__1 = k + (k - 1) * a_dim1;
#line 369 "zsytri.f"
		i__2 = kp + (k - 1) * a_dim1;
#line 369 "zsytri.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 370 "zsytri.f"
		i__1 = kp + (k - 1) * a_dim1;
#line 370 "zsytri.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 371 "zsytri.f"
	    }
#line 372 "zsytri.f"
	}

#line 374 "zsytri.f"
	k -= kstep;
#line 375 "zsytri.f"
	goto L50;
#line 376 "zsytri.f"
L60:
#line 377 "zsytri.f"
	;
#line 377 "zsytri.f"
    }

#line 379 "zsytri.f"
    return 0;

/*     End of ZSYTRI */

} /* zsytri_ */

