#line 1 "zhptri.f"
/* zhptri.f -- translated by f2c (version 20100827).
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

#line 1 "zhptri.f"
/* Table of constant values */

static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b ZHPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPTRI( UPLO, N, AP, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPTRI computes the inverse of a complex Hermitian indefinite matrix */
/* > A in packed storage using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by ZHPTRF. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the block diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by ZHPTRF, */
/* >          stored as a packed triangular matrix. */
/* > */
/* >          On exit, if INFO = 0, the (Hermitian) inverse of the original */
/* >          matrix, stored as a packed triangular matrix. The j-th column */
/* >          of inv(A) is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', */
/* >             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZHPTRF. */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublecomplex *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer kc, kp, kx, kpc, npp;
    static doublereal akp1;
    static doublecomplex temp, akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zhpmv_(char *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), zswap_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
    static integer kcnext;


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

#line 155 "zhptri.f"
    /* Parameter adjustments */
#line 155 "zhptri.f"
    --work;
#line 155 "zhptri.f"
    --ipiv;
#line 155 "zhptri.f"
    --ap;
#line 155 "zhptri.f"

#line 155 "zhptri.f"
    /* Function Body */
#line 155 "zhptri.f"
    *info = 0;
#line 156 "zhptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 157 "zhptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 158 "zhptri.f"
	*info = -1;
#line 159 "zhptri.f"
    } else if (*n < 0) {
#line 160 "zhptri.f"
	*info = -2;
#line 161 "zhptri.f"
    }
#line 162 "zhptri.f"
    if (*info != 0) {
#line 163 "zhptri.f"
	i__1 = -(*info);
#line 163 "zhptri.f"
	xerbla_("ZHPTRI", &i__1, (ftnlen)6);
#line 164 "zhptri.f"
	return 0;
#line 165 "zhptri.f"
    }

/*     Quick return if possible */

#line 169 "zhptri.f"
    if (*n == 0) {
#line 169 "zhptri.f"
	return 0;
#line 169 "zhptri.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 174 "zhptri.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 178 "zhptri.f"
	kp = *n * (*n + 1) / 2;
#line 179 "zhptri.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 180 "zhptri.f"
	    i__1 = kp;
#line 180 "zhptri.f"
	    if (ipiv[*info] > 0 && (ap[i__1].r == 0. && ap[i__1].i == 0.)) {
#line 180 "zhptri.f"
		return 0;
#line 180 "zhptri.f"
	    }
#line 182 "zhptri.f"
	    kp -= *info;
#line 183 "zhptri.f"
/* L10: */
#line 183 "zhptri.f"
	}
#line 184 "zhptri.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 188 "zhptri.f"
	kp = 1;
#line 189 "zhptri.f"
	i__1 = *n;
#line 189 "zhptri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 190 "zhptri.f"
	    i__2 = kp;
#line 190 "zhptri.f"
	    if (ipiv[*info] > 0 && (ap[i__2].r == 0. && ap[i__2].i == 0.)) {
#line 190 "zhptri.f"
		return 0;
#line 190 "zhptri.f"
	    }
#line 192 "zhptri.f"
	    kp = kp + *n - *info + 1;
#line 193 "zhptri.f"
/* L20: */
#line 193 "zhptri.f"
	}
#line 194 "zhptri.f"
    }
#line 195 "zhptri.f"
    *info = 0;

#line 197 "zhptri.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**H. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 204 "zhptri.f"
	k = 1;
#line 205 "zhptri.f"
	kc = 1;
#line 206 "zhptri.f"
L30:

/*        If K > N, exit from loop. */

#line 210 "zhptri.f"
	if (k > *n) {
#line 210 "zhptri.f"
	    goto L50;
#line 210 "zhptri.f"
	}

#line 213 "zhptri.f"
	kcnext = kc + k;
#line 214 "zhptri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 220 "zhptri.f"
	    i__1 = kc + k - 1;
#line 220 "zhptri.f"
	    i__2 = kc + k - 1;
#line 220 "zhptri.f"
	    d__1 = 1. / ap[i__2].r;
#line 220 "zhptri.f"
	    ap[i__1].r = d__1, ap[i__1].i = 0.;

/*           Compute column K of the inverse. */

#line 224 "zhptri.f"
	    if (k > 1) {
#line 225 "zhptri.f"
		i__1 = k - 1;
#line 225 "zhptri.f"
		zcopy_(&i__1, &ap[kc], &c__1, &work[1], &c__1);
#line 226 "zhptri.f"
		i__1 = k - 1;
#line 226 "zhptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 226 "zhptri.f"
		zhpmv_(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &
			ap[kc], &c__1, (ftnlen)1);
#line 228 "zhptri.f"
		i__1 = kc + k - 1;
#line 228 "zhptri.f"
		i__2 = kc + k - 1;
#line 228 "zhptri.f"
		i__3 = k - 1;
#line 228 "zhptri.f"
		zdotc_(&z__2, &i__3, &work[1], &c__1, &ap[kc], &c__1);
#line 228 "zhptri.f"
		d__1 = z__2.r;
#line 228 "zhptri.f"
		z__1.r = ap[i__2].r - d__1, z__1.i = ap[i__2].i;
#line 228 "zhptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 230 "zhptri.f"
	    }
#line 231 "zhptri.f"
	    kstep = 1;
#line 232 "zhptri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 238 "zhptri.f"
	    t = z_abs(&ap[kcnext + k - 1]);
#line 239 "zhptri.f"
	    i__1 = kc + k - 1;
#line 239 "zhptri.f"
	    ak = ap[i__1].r / t;
#line 240 "zhptri.f"
	    i__1 = kcnext + k;
#line 240 "zhptri.f"
	    akp1 = ap[i__1].r / t;
#line 241 "zhptri.f"
	    i__1 = kcnext + k - 1;
#line 241 "zhptri.f"
	    z__1.r = ap[i__1].r / t, z__1.i = ap[i__1].i / t;
#line 241 "zhptri.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 242 "zhptri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 243 "zhptri.f"
	    i__1 = kc + k - 1;
#line 243 "zhptri.f"
	    d__1 = akp1 / d__;
#line 243 "zhptri.f"
	    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 244 "zhptri.f"
	    i__1 = kcnext + k;
#line 244 "zhptri.f"
	    d__1 = ak / d__;
#line 244 "zhptri.f"
	    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 245 "zhptri.f"
	    i__1 = kcnext + k - 1;
#line 245 "zhptri.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 245 "zhptri.f"
	    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
#line 245 "zhptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;

/*           Compute columns K and K+1 of the inverse. */

#line 249 "zhptri.f"
	    if (k > 1) {
#line 250 "zhptri.f"
		i__1 = k - 1;
#line 250 "zhptri.f"
		zcopy_(&i__1, &ap[kc], &c__1, &work[1], &c__1);
#line 251 "zhptri.f"
		i__1 = k - 1;
#line 251 "zhptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 251 "zhptri.f"
		zhpmv_(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &
			ap[kc], &c__1, (ftnlen)1);
#line 253 "zhptri.f"
		i__1 = kc + k - 1;
#line 253 "zhptri.f"
		i__2 = kc + k - 1;
#line 253 "zhptri.f"
		i__3 = k - 1;
#line 253 "zhptri.f"
		zdotc_(&z__2, &i__3, &work[1], &c__1, &ap[kc], &c__1);
#line 253 "zhptri.f"
		d__1 = z__2.r;
#line 253 "zhptri.f"
		z__1.r = ap[i__2].r - d__1, z__1.i = ap[i__2].i;
#line 253 "zhptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 255 "zhptri.f"
		i__1 = kcnext + k - 1;
#line 255 "zhptri.f"
		i__2 = kcnext + k - 1;
#line 255 "zhptri.f"
		i__3 = k - 1;
#line 255 "zhptri.f"
		zdotc_(&z__2, &i__3, &ap[kc], &c__1, &ap[kcnext], &c__1);
#line 255 "zhptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 255 "zhptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 258 "zhptri.f"
		i__1 = k - 1;
#line 258 "zhptri.f"
		zcopy_(&i__1, &ap[kcnext], &c__1, &work[1], &c__1);
#line 259 "zhptri.f"
		i__1 = k - 1;
#line 259 "zhptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 259 "zhptri.f"
		zhpmv_(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &
			ap[kcnext], &c__1, (ftnlen)1);
#line 261 "zhptri.f"
		i__1 = kcnext + k;
#line 261 "zhptri.f"
		i__2 = kcnext + k;
#line 261 "zhptri.f"
		i__3 = k - 1;
#line 261 "zhptri.f"
		zdotc_(&z__2, &i__3, &work[1], &c__1, &ap[kcnext], &c__1);
#line 261 "zhptri.f"
		d__1 = z__2.r;
#line 261 "zhptri.f"
		z__1.r = ap[i__2].r - d__1, z__1.i = ap[i__2].i;
#line 261 "zhptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 264 "zhptri.f"
	    }
#line 265 "zhptri.f"
	    kstep = 2;
#line 266 "zhptri.f"
	    kcnext = kcnext + k + 1;
#line 267 "zhptri.f"
	}

#line 269 "zhptri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 270 "zhptri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 275 "zhptri.f"
	    kpc = (kp - 1) * kp / 2 + 1;
#line 276 "zhptri.f"
	    i__1 = kp - 1;
#line 276 "zhptri.f"
	    zswap_(&i__1, &ap[kc], &c__1, &ap[kpc], &c__1);
#line 277 "zhptri.f"
	    kx = kpc + kp - 1;
#line 278 "zhptri.f"
	    i__1 = k - 1;
#line 278 "zhptri.f"
	    for (j = kp + 1; j <= i__1; ++j) {
#line 279 "zhptri.f"
		kx = kx + j - 1;
#line 280 "zhptri.f"
		d_cnjg(&z__1, &ap[kc + j - 1]);
#line 280 "zhptri.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 281 "zhptri.f"
		i__2 = kc + j - 1;
#line 281 "zhptri.f"
		d_cnjg(&z__1, &ap[kx]);
#line 281 "zhptri.f"
		ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 282 "zhptri.f"
		i__2 = kx;
#line 282 "zhptri.f"
		ap[i__2].r = temp.r, ap[i__2].i = temp.i;
#line 283 "zhptri.f"
/* L40: */
#line 283 "zhptri.f"
	    }
#line 284 "zhptri.f"
	    i__1 = kc + kp - 1;
#line 284 "zhptri.f"
	    d_cnjg(&z__1, &ap[kc + kp - 1]);
#line 284 "zhptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 285 "zhptri.f"
	    i__1 = kc + k - 1;
#line 285 "zhptri.f"
	    temp.r = ap[i__1].r, temp.i = ap[i__1].i;
#line 286 "zhptri.f"
	    i__1 = kc + k - 1;
#line 286 "zhptri.f"
	    i__2 = kpc + kp - 1;
#line 286 "zhptri.f"
	    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 287 "zhptri.f"
	    i__1 = kpc + kp - 1;
#line 287 "zhptri.f"
	    ap[i__1].r = temp.r, ap[i__1].i = temp.i;
#line 288 "zhptri.f"
	    if (kstep == 2) {
#line 289 "zhptri.f"
		i__1 = kc + k + k - 1;
#line 289 "zhptri.f"
		temp.r = ap[i__1].r, temp.i = ap[i__1].i;
#line 290 "zhptri.f"
		i__1 = kc + k + k - 1;
#line 290 "zhptri.f"
		i__2 = kc + k + kp - 1;
#line 290 "zhptri.f"
		ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 291 "zhptri.f"
		i__1 = kc + k + kp - 1;
#line 291 "zhptri.f"
		ap[i__1].r = temp.r, ap[i__1].i = temp.i;
#line 292 "zhptri.f"
	    }
#line 293 "zhptri.f"
	}

#line 295 "zhptri.f"
	k += kstep;
#line 296 "zhptri.f"
	kc = kcnext;
#line 297 "zhptri.f"
	goto L30;
#line 298 "zhptri.f"
L50:

#line 300 "zhptri.f"
	;
#line 300 "zhptri.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**H. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 307 "zhptri.f"
	npp = *n * (*n + 1) / 2;
#line 308 "zhptri.f"
	k = *n;
#line 309 "zhptri.f"
	kc = npp;
#line 310 "zhptri.f"
L60:

/*        If K < 1, exit from loop. */

#line 314 "zhptri.f"
	if (k < 1) {
#line 314 "zhptri.f"
	    goto L80;
#line 314 "zhptri.f"
	}

#line 317 "zhptri.f"
	kcnext = kc - (*n - k + 2);
#line 318 "zhptri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 324 "zhptri.f"
	    i__1 = kc;
#line 324 "zhptri.f"
	    i__2 = kc;
#line 324 "zhptri.f"
	    d__1 = 1. / ap[i__2].r;
#line 324 "zhptri.f"
	    ap[i__1].r = d__1, ap[i__1].i = 0.;

/*           Compute column K of the inverse. */

#line 328 "zhptri.f"
	    if (k < *n) {
#line 329 "zhptri.f"
		i__1 = *n - k;
#line 329 "zhptri.f"
		zcopy_(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
#line 330 "zhptri.f"
		i__1 = *n - k;
#line 330 "zhptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 330 "zhptri.f"
		zhpmv_(uplo, &i__1, &z__1, &ap[kc + *n - k + 1], &work[1], &
			c__1, &c_b2, &ap[kc + 1], &c__1, (ftnlen)1);
#line 332 "zhptri.f"
		i__1 = kc;
#line 332 "zhptri.f"
		i__2 = kc;
#line 332 "zhptri.f"
		i__3 = *n - k;
#line 332 "zhptri.f"
		zdotc_(&z__2, &i__3, &work[1], &c__1, &ap[kc + 1], &c__1);
#line 332 "zhptri.f"
		d__1 = z__2.r;
#line 332 "zhptri.f"
		z__1.r = ap[i__2].r - d__1, z__1.i = ap[i__2].i;
#line 332 "zhptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 334 "zhptri.f"
	    }
#line 335 "zhptri.f"
	    kstep = 1;
#line 336 "zhptri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 342 "zhptri.f"
	    t = z_abs(&ap[kcnext + 1]);
#line 343 "zhptri.f"
	    i__1 = kcnext;
#line 343 "zhptri.f"
	    ak = ap[i__1].r / t;
#line 344 "zhptri.f"
	    i__1 = kc;
#line 344 "zhptri.f"
	    akp1 = ap[i__1].r / t;
#line 345 "zhptri.f"
	    i__1 = kcnext + 1;
#line 345 "zhptri.f"
	    z__1.r = ap[i__1].r / t, z__1.i = ap[i__1].i / t;
#line 345 "zhptri.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 346 "zhptri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 347 "zhptri.f"
	    i__1 = kcnext;
#line 347 "zhptri.f"
	    d__1 = akp1 / d__;
#line 347 "zhptri.f"
	    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 348 "zhptri.f"
	    i__1 = kc;
#line 348 "zhptri.f"
	    d__1 = ak / d__;
#line 348 "zhptri.f"
	    ap[i__1].r = d__1, ap[i__1].i = 0.;
#line 349 "zhptri.f"
	    i__1 = kcnext + 1;
#line 349 "zhptri.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 349 "zhptri.f"
	    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
#line 349 "zhptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;

/*           Compute columns K-1 and K of the inverse. */

#line 353 "zhptri.f"
	    if (k < *n) {
#line 354 "zhptri.f"
		i__1 = *n - k;
#line 354 "zhptri.f"
		zcopy_(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
#line 355 "zhptri.f"
		i__1 = *n - k;
#line 355 "zhptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 355 "zhptri.f"
		zhpmv_(uplo, &i__1, &z__1, &ap[kc + (*n - k + 1)], &work[1], &
			c__1, &c_b2, &ap[kc + 1], &c__1, (ftnlen)1);
#line 357 "zhptri.f"
		i__1 = kc;
#line 357 "zhptri.f"
		i__2 = kc;
#line 357 "zhptri.f"
		i__3 = *n - k;
#line 357 "zhptri.f"
		zdotc_(&z__2, &i__3, &work[1], &c__1, &ap[kc + 1], &c__1);
#line 357 "zhptri.f"
		d__1 = z__2.r;
#line 357 "zhptri.f"
		z__1.r = ap[i__2].r - d__1, z__1.i = ap[i__2].i;
#line 357 "zhptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 359 "zhptri.f"
		i__1 = kcnext + 1;
#line 359 "zhptri.f"
		i__2 = kcnext + 1;
#line 359 "zhptri.f"
		i__3 = *n - k;
#line 359 "zhptri.f"
		zdotc_(&z__2, &i__3, &ap[kc + 1], &c__1, &ap[kcnext + 2], &
			c__1);
#line 359 "zhptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 359 "zhptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 362 "zhptri.f"
		i__1 = *n - k;
#line 362 "zhptri.f"
		zcopy_(&i__1, &ap[kcnext + 2], &c__1, &work[1], &c__1);
#line 363 "zhptri.f"
		i__1 = *n - k;
#line 363 "zhptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 363 "zhptri.f"
		zhpmv_(uplo, &i__1, &z__1, &ap[kc + (*n - k + 1)], &work[1], &
			c__1, &c_b2, &ap[kcnext + 2], &c__1, (ftnlen)1);
#line 365 "zhptri.f"
		i__1 = kcnext;
#line 365 "zhptri.f"
		i__2 = kcnext;
#line 365 "zhptri.f"
		i__3 = *n - k;
#line 365 "zhptri.f"
		zdotc_(&z__2, &i__3, &work[1], &c__1, &ap[kcnext + 2], &c__1);
#line 365 "zhptri.f"
		d__1 = z__2.r;
#line 365 "zhptri.f"
		z__1.r = ap[i__2].r - d__1, z__1.i = ap[i__2].i;
#line 365 "zhptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 368 "zhptri.f"
	    }
#line 369 "zhptri.f"
	    kstep = 2;
#line 370 "zhptri.f"
	    kcnext -= *n - k + 3;
#line 371 "zhptri.f"
	}

#line 373 "zhptri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 374 "zhptri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 379 "zhptri.f"
	    kpc = npp - (*n - kp + 1) * (*n - kp + 2) / 2 + 1;
#line 380 "zhptri.f"
	    if (kp < *n) {
#line 380 "zhptri.f"
		i__1 = *n - kp;
#line 380 "zhptri.f"
		zswap_(&i__1, &ap[kc + kp - k + 1], &c__1, &ap[kpc + 1], &
			c__1);
#line 380 "zhptri.f"
	    }
#line 382 "zhptri.f"
	    kx = kc + kp - k;
#line 383 "zhptri.f"
	    i__1 = kp - 1;
#line 383 "zhptri.f"
	    for (j = k + 1; j <= i__1; ++j) {
#line 384 "zhptri.f"
		kx = kx + *n - j + 1;
#line 385 "zhptri.f"
		d_cnjg(&z__1, &ap[kc + j - k]);
#line 385 "zhptri.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 386 "zhptri.f"
		i__2 = kc + j - k;
#line 386 "zhptri.f"
		d_cnjg(&z__1, &ap[kx]);
#line 386 "zhptri.f"
		ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 387 "zhptri.f"
		i__2 = kx;
#line 387 "zhptri.f"
		ap[i__2].r = temp.r, ap[i__2].i = temp.i;
#line 388 "zhptri.f"
/* L70: */
#line 388 "zhptri.f"
	    }
#line 389 "zhptri.f"
	    i__1 = kc + kp - k;
#line 389 "zhptri.f"
	    d_cnjg(&z__1, &ap[kc + kp - k]);
#line 389 "zhptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 390 "zhptri.f"
	    i__1 = kc;
#line 390 "zhptri.f"
	    temp.r = ap[i__1].r, temp.i = ap[i__1].i;
#line 391 "zhptri.f"
	    i__1 = kc;
#line 391 "zhptri.f"
	    i__2 = kpc;
#line 391 "zhptri.f"
	    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 392 "zhptri.f"
	    i__1 = kpc;
#line 392 "zhptri.f"
	    ap[i__1].r = temp.r, ap[i__1].i = temp.i;
#line 393 "zhptri.f"
	    if (kstep == 2) {
#line 394 "zhptri.f"
		i__1 = kc - *n + k - 1;
#line 394 "zhptri.f"
		temp.r = ap[i__1].r, temp.i = ap[i__1].i;
#line 395 "zhptri.f"
		i__1 = kc - *n + k - 1;
#line 395 "zhptri.f"
		i__2 = kc - *n + kp - 1;
#line 395 "zhptri.f"
		ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 396 "zhptri.f"
		i__1 = kc - *n + kp - 1;
#line 396 "zhptri.f"
		ap[i__1].r = temp.r, ap[i__1].i = temp.i;
#line 397 "zhptri.f"
	    }
#line 398 "zhptri.f"
	}

#line 400 "zhptri.f"
	k -= kstep;
#line 401 "zhptri.f"
	kc = kcnext;
#line 402 "zhptri.f"
	goto L60;
#line 403 "zhptri.f"
L80:
#line 404 "zhptri.f"
	;
#line 404 "zhptri.f"
    }

#line 406 "zhptri.f"
    return 0;

/*     End of ZHPTRI */

} /* zhptri_ */

