#line 1 "chetri.f"
/* chetri.f -- translated by f2c (version 20100827).
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

#line 1 "chetri.f"
/* Table of constant values */

static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b CHETRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHETRI computes the inverse of a complex Hermitian indefinite matrix */
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
/* >          On entry, the block diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by CHETRF. */
/* > */
/* >          On exit, if INFO = 0, the (Hermitian) inverse of the original */
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
/* >          as determined by CHETRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N) */
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
/* Subroutine */ int chetri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer kp;
    static doublereal akp1;
    static doublecomplex temp, akkp1;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int chemv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), ccopy_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    , cswap_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

#line 160 "chetri.f"
    /* Parameter adjustments */
#line 160 "chetri.f"
    a_dim1 = *lda;
#line 160 "chetri.f"
    a_offset = 1 + a_dim1;
#line 160 "chetri.f"
    a -= a_offset;
#line 160 "chetri.f"
    --ipiv;
#line 160 "chetri.f"
    --work;
#line 160 "chetri.f"

#line 160 "chetri.f"
    /* Function Body */
#line 160 "chetri.f"
    *info = 0;
#line 161 "chetri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 162 "chetri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 163 "chetri.f"
	*info = -1;
#line 164 "chetri.f"
    } else if (*n < 0) {
#line 165 "chetri.f"
	*info = -2;
#line 166 "chetri.f"
    } else if (*lda < max(1,*n)) {
#line 167 "chetri.f"
	*info = -4;
#line 168 "chetri.f"
    }
#line 169 "chetri.f"
    if (*info != 0) {
#line 170 "chetri.f"
	i__1 = -(*info);
#line 170 "chetri.f"
	xerbla_("CHETRI", &i__1, (ftnlen)6);
#line 171 "chetri.f"
	return 0;
#line 172 "chetri.f"
    }

/*     Quick return if possible */

#line 176 "chetri.f"
    if (*n == 0) {
#line 176 "chetri.f"
	return 0;
#line 176 "chetri.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 181 "chetri.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 185 "chetri.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 186 "chetri.f"
	    i__1 = *info + *info * a_dim1;
#line 186 "chetri.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 186 "chetri.f"
		return 0;
#line 186 "chetri.f"
	    }
#line 188 "chetri.f"
/* L10: */
#line 188 "chetri.f"
	}
#line 189 "chetri.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 193 "chetri.f"
	i__1 = *n;
#line 193 "chetri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 194 "chetri.f"
	    i__2 = *info + *info * a_dim1;
#line 194 "chetri.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 194 "chetri.f"
		return 0;
#line 194 "chetri.f"
	    }
#line 196 "chetri.f"
/* L20: */
#line 196 "chetri.f"
	}
#line 197 "chetri.f"
    }
#line 198 "chetri.f"
    *info = 0;

#line 200 "chetri.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**H. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 207 "chetri.f"
	k = 1;
#line 208 "chetri.f"
L30:

/*        If K > N, exit from loop. */

#line 212 "chetri.f"
	if (k > *n) {
#line 212 "chetri.f"
	    goto L50;
#line 212 "chetri.f"
	}

#line 215 "chetri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 221 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 221 "chetri.f"
	    i__2 = k + k * a_dim1;
#line 221 "chetri.f"
	    d__1 = 1. / a[i__2].r;
#line 221 "chetri.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;

/*           Compute column K of the inverse. */

#line 225 "chetri.f"
	    if (k > 1) {
#line 226 "chetri.f"
		i__1 = k - 1;
#line 226 "chetri.f"
		ccopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 227 "chetri.f"
		i__1 = k - 1;
#line 227 "chetri.f"
		z__1.r = -1., z__1.i = -0.;
#line 227 "chetri.f"
		chemv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 229 "chetri.f"
		i__1 = k + k * a_dim1;
#line 229 "chetri.f"
		i__2 = k + k * a_dim1;
#line 229 "chetri.f"
		i__3 = k - 1;
#line 229 "chetri.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 229 "chetri.f"
		d__1 = z__2.r;
#line 229 "chetri.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 229 "chetri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 231 "chetri.f"
	    }
#line 232 "chetri.f"
	    kstep = 1;
#line 233 "chetri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 239 "chetri.f"
	    t = z_abs(&a[k + (k + 1) * a_dim1]);
#line 240 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 240 "chetri.f"
	    ak = a[i__1].r / t;
#line 241 "chetri.f"
	    i__1 = k + 1 + (k + 1) * a_dim1;
#line 241 "chetri.f"
	    akp1 = a[i__1].r / t;
#line 242 "chetri.f"
	    i__1 = k + (k + 1) * a_dim1;
#line 242 "chetri.f"
	    z__1.r = a[i__1].r / t, z__1.i = a[i__1].i / t;
#line 242 "chetri.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 243 "chetri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 244 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 244 "chetri.f"
	    d__1 = akp1 / d__;
#line 244 "chetri.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 245 "chetri.f"
	    i__1 = k + 1 + (k + 1) * a_dim1;
#line 245 "chetri.f"
	    d__1 = ak / d__;
#line 245 "chetri.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 246 "chetri.f"
	    i__1 = k + (k + 1) * a_dim1;
#line 246 "chetri.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 246 "chetri.f"
	    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
#line 246 "chetri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute columns K and K+1 of the inverse. */

#line 250 "chetri.f"
	    if (k > 1) {
#line 251 "chetri.f"
		i__1 = k - 1;
#line 251 "chetri.f"
		ccopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 252 "chetri.f"
		i__1 = k - 1;
#line 252 "chetri.f"
		z__1.r = -1., z__1.i = -0.;
#line 252 "chetri.f"
		chemv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 254 "chetri.f"
		i__1 = k + k * a_dim1;
#line 254 "chetri.f"
		i__2 = k + k * a_dim1;
#line 254 "chetri.f"
		i__3 = k - 1;
#line 254 "chetri.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 254 "chetri.f"
		d__1 = z__2.r;
#line 254 "chetri.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 254 "chetri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 256 "chetri.f"
		i__1 = k + (k + 1) * a_dim1;
#line 256 "chetri.f"
		i__2 = k + (k + 1) * a_dim1;
#line 256 "chetri.f"
		i__3 = k - 1;
#line 256 "chetri.f"
		cdotc_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * 
			a_dim1 + 1], &c__1);
#line 256 "chetri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 256 "chetri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 258 "chetri.f"
		i__1 = k - 1;
#line 258 "chetri.f"
		ccopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &
			c__1);
#line 259 "chetri.f"
		i__1 = k - 1;
#line 259 "chetri.f"
		z__1.r = -1., z__1.i = -0.;
#line 259 "chetri.f"
		chemv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[(k + 1) * a_dim1 + 1], &c__1, (ftnlen)1);
#line 261 "chetri.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 261 "chetri.f"
		i__2 = k + 1 + (k + 1) * a_dim1;
#line 261 "chetri.f"
		i__3 = k - 1;
#line 261 "chetri.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1]
			, &c__1);
#line 261 "chetri.f"
		d__1 = z__2.r;
#line 261 "chetri.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 261 "chetri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 264 "chetri.f"
	    }
#line 265 "chetri.f"
	    kstep = 2;
#line 266 "chetri.f"
	}

#line 268 "chetri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 269 "chetri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 274 "chetri.f"
	    i__1 = kp - 1;
#line 274 "chetri.f"
	    cswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
		    c__1);
#line 275 "chetri.f"
	    i__1 = k - 1;
#line 275 "chetri.f"
	    for (j = kp + 1; j <= i__1; ++j) {
#line 276 "chetri.f"
		d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 276 "chetri.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 277 "chetri.f"
		i__2 = j + k * a_dim1;
#line 277 "chetri.f"
		d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 277 "chetri.f"
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 278 "chetri.f"
		i__2 = kp + j * a_dim1;
#line 278 "chetri.f"
		a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 279 "chetri.f"
/* L40: */
#line 279 "chetri.f"
	    }
#line 280 "chetri.f"
	    i__1 = kp + k * a_dim1;
#line 280 "chetri.f"
	    d_cnjg(&z__1, &a[kp + k * a_dim1]);
#line 280 "chetri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 281 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 281 "chetri.f"
	    temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 282 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 282 "chetri.f"
	    i__2 = kp + kp * a_dim1;
#line 282 "chetri.f"
	    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 283 "chetri.f"
	    i__1 = kp + kp * a_dim1;
#line 283 "chetri.f"
	    a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 284 "chetri.f"
	    if (kstep == 2) {
#line 285 "chetri.f"
		i__1 = k + (k + 1) * a_dim1;
#line 285 "chetri.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 286 "chetri.f"
		i__1 = k + (k + 1) * a_dim1;
#line 286 "chetri.f"
		i__2 = kp + (k + 1) * a_dim1;
#line 286 "chetri.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 287 "chetri.f"
		i__1 = kp + (k + 1) * a_dim1;
#line 287 "chetri.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 288 "chetri.f"
	    }
#line 289 "chetri.f"
	}

#line 291 "chetri.f"
	k += kstep;
#line 292 "chetri.f"
	goto L30;
#line 293 "chetri.f"
L50:

#line 295 "chetri.f"
	;
#line 295 "chetri.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**H. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 302 "chetri.f"
	k = *n;
#line 303 "chetri.f"
L60:

/*        If K < 1, exit from loop. */

#line 307 "chetri.f"
	if (k < 1) {
#line 307 "chetri.f"
	    goto L80;
#line 307 "chetri.f"
	}

#line 310 "chetri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 316 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 316 "chetri.f"
	    i__2 = k + k * a_dim1;
#line 316 "chetri.f"
	    d__1 = 1. / a[i__2].r;
#line 316 "chetri.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;

/*           Compute column K of the inverse. */

#line 320 "chetri.f"
	    if (k < *n) {
#line 321 "chetri.f"
		i__1 = *n - k;
#line 321 "chetri.f"
		ccopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 322 "chetri.f"
		i__1 = *n - k;
#line 322 "chetri.f"
		z__1.r = -1., z__1.i = -0.;
#line 322 "chetri.f"
		chemv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1,
			 (ftnlen)1);
#line 324 "chetri.f"
		i__1 = k + k * a_dim1;
#line 324 "chetri.f"
		i__2 = k + k * a_dim1;
#line 324 "chetri.f"
		i__3 = *n - k;
#line 324 "chetri.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], 
			&c__1);
#line 324 "chetri.f"
		d__1 = z__2.r;
#line 324 "chetri.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 324 "chetri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 326 "chetri.f"
	    }
#line 327 "chetri.f"
	    kstep = 1;
#line 328 "chetri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 334 "chetri.f"
	    t = z_abs(&a[k + (k - 1) * a_dim1]);
#line 335 "chetri.f"
	    i__1 = k - 1 + (k - 1) * a_dim1;
#line 335 "chetri.f"
	    ak = a[i__1].r / t;
#line 336 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 336 "chetri.f"
	    akp1 = a[i__1].r / t;
#line 337 "chetri.f"
	    i__1 = k + (k - 1) * a_dim1;
#line 337 "chetri.f"
	    z__1.r = a[i__1].r / t, z__1.i = a[i__1].i / t;
#line 337 "chetri.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 338 "chetri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 339 "chetri.f"
	    i__1 = k - 1 + (k - 1) * a_dim1;
#line 339 "chetri.f"
	    d__1 = akp1 / d__;
#line 339 "chetri.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 340 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 340 "chetri.f"
	    d__1 = ak / d__;
#line 340 "chetri.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 341 "chetri.f"
	    i__1 = k + (k - 1) * a_dim1;
#line 341 "chetri.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 341 "chetri.f"
	    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
#line 341 "chetri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute columns K-1 and K of the inverse. */

#line 345 "chetri.f"
	    if (k < *n) {
#line 346 "chetri.f"
		i__1 = *n - k;
#line 346 "chetri.f"
		ccopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 347 "chetri.f"
		i__1 = *n - k;
#line 347 "chetri.f"
		z__1.r = -1., z__1.i = -0.;
#line 347 "chetri.f"
		chemv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1,
			 (ftnlen)1);
#line 349 "chetri.f"
		i__1 = k + k * a_dim1;
#line 349 "chetri.f"
		i__2 = k + k * a_dim1;
#line 349 "chetri.f"
		i__3 = *n - k;
#line 349 "chetri.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], 
			&c__1);
#line 349 "chetri.f"
		d__1 = z__2.r;
#line 349 "chetri.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 349 "chetri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 351 "chetri.f"
		i__1 = k + (k - 1) * a_dim1;
#line 351 "chetri.f"
		i__2 = k + (k - 1) * a_dim1;
#line 351 "chetri.f"
		i__3 = *n - k;
#line 351 "chetri.f"
		cdotc_(&z__2, &i__3, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 
			+ (k - 1) * a_dim1], &c__1);
#line 351 "chetri.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 351 "chetri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 354 "chetri.f"
		i__1 = *n - k;
#line 354 "chetri.f"
		ccopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &
			c__1);
#line 355 "chetri.f"
		i__1 = *n - k;
#line 355 "chetri.f"
		z__1.r = -1., z__1.i = -0.;
#line 355 "chetri.f"
		chemv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + (k - 1) * a_dim1], 
			&c__1, (ftnlen)1);
#line 357 "chetri.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 357 "chetri.f"
		i__2 = k - 1 + (k - 1) * a_dim1;
#line 357 "chetri.f"
		i__3 = *n - k;
#line 357 "chetri.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + (k - 1) * 
			a_dim1], &c__1);
#line 357 "chetri.f"
		d__1 = z__2.r;
#line 357 "chetri.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 357 "chetri.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 360 "chetri.f"
	    }
#line 361 "chetri.f"
	    kstep = 2;
#line 362 "chetri.f"
	}

#line 364 "chetri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 365 "chetri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 370 "chetri.f"
	    if (kp < *n) {
#line 370 "chetri.f"
		i__1 = *n - kp;
#line 370 "chetri.f"
		cswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp *
			 a_dim1], &c__1);
#line 370 "chetri.f"
	    }
#line 372 "chetri.f"
	    i__1 = kp - 1;
#line 372 "chetri.f"
	    for (j = k + 1; j <= i__1; ++j) {
#line 373 "chetri.f"
		d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 373 "chetri.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 374 "chetri.f"
		i__2 = j + k * a_dim1;
#line 374 "chetri.f"
		d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 374 "chetri.f"
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 375 "chetri.f"
		i__2 = kp + j * a_dim1;
#line 375 "chetri.f"
		a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 376 "chetri.f"
/* L70: */
#line 376 "chetri.f"
	    }
#line 377 "chetri.f"
	    i__1 = kp + k * a_dim1;
#line 377 "chetri.f"
	    d_cnjg(&z__1, &a[kp + k * a_dim1]);
#line 377 "chetri.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 378 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 378 "chetri.f"
	    temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 379 "chetri.f"
	    i__1 = k + k * a_dim1;
#line 379 "chetri.f"
	    i__2 = kp + kp * a_dim1;
#line 379 "chetri.f"
	    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 380 "chetri.f"
	    i__1 = kp + kp * a_dim1;
#line 380 "chetri.f"
	    a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 381 "chetri.f"
	    if (kstep == 2) {
#line 382 "chetri.f"
		i__1 = k + (k - 1) * a_dim1;
#line 382 "chetri.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 383 "chetri.f"
		i__1 = k + (k - 1) * a_dim1;
#line 383 "chetri.f"
		i__2 = kp + (k - 1) * a_dim1;
#line 383 "chetri.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 384 "chetri.f"
		i__1 = kp + (k - 1) * a_dim1;
#line 384 "chetri.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 385 "chetri.f"
	    }
#line 386 "chetri.f"
	}

#line 388 "chetri.f"
	k -= kstep;
#line 389 "chetri.f"
	goto L60;
#line 390 "chetri.f"
L80:
#line 391 "chetri.f"
	;
#line 391 "chetri.f"
    }

#line 393 "chetri.f"
    return 0;

/*     End of CHETRI */

} /* chetri_ */

