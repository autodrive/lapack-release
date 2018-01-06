#line 1 "csptri.f"
/* csptri.f -- translated by f2c (version 20100827).
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

#line 1 "csptri.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b CSPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSPTRI( UPLO, N, AP, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSPTRI computes the inverse of a complex symmetric indefinite matrix */
/* > A in packed storage using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by CSPTRF. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the block diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by CSPTRF, */
/* >          stored as a packed triangular matrix. */
/* > */
/* >          On exit, if INFO = 0, the (symmetric) inverse of the original */
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
/* >          as determined by CSPTRF. */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int csptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublecomplex *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex d__;
    static integer j, k;
    static doublecomplex t, ak;
    static integer kc, kp, kx, kpc, npp;
    static doublecomplex akp1, temp, akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer kstep;
    extern /* Subroutine */ int cspmv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer kcnext;


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

#line 153 "csptri.f"
    /* Parameter adjustments */
#line 153 "csptri.f"
    --work;
#line 153 "csptri.f"
    --ipiv;
#line 153 "csptri.f"
    --ap;
#line 153 "csptri.f"

#line 153 "csptri.f"
    /* Function Body */
#line 153 "csptri.f"
    *info = 0;
#line 154 "csptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 155 "csptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 156 "csptri.f"
	*info = -1;
#line 157 "csptri.f"
    } else if (*n < 0) {
#line 158 "csptri.f"
	*info = -2;
#line 159 "csptri.f"
    }
#line 160 "csptri.f"
    if (*info != 0) {
#line 161 "csptri.f"
	i__1 = -(*info);
#line 161 "csptri.f"
	xerbla_("CSPTRI", &i__1, (ftnlen)6);
#line 162 "csptri.f"
	return 0;
#line 163 "csptri.f"
    }

/*     Quick return if possible */

#line 167 "csptri.f"
    if (*n == 0) {
#line 167 "csptri.f"
	return 0;
#line 167 "csptri.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 172 "csptri.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 176 "csptri.f"
	kp = *n * (*n + 1) / 2;
#line 177 "csptri.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 178 "csptri.f"
	    i__1 = kp;
#line 178 "csptri.f"
	    if (ipiv[*info] > 0 && (ap[i__1].r == 0. && ap[i__1].i == 0.)) {
#line 178 "csptri.f"
		return 0;
#line 178 "csptri.f"
	    }
#line 180 "csptri.f"
	    kp -= *info;
#line 181 "csptri.f"
/* L10: */
#line 181 "csptri.f"
	}
#line 182 "csptri.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 186 "csptri.f"
	kp = 1;
#line 187 "csptri.f"
	i__1 = *n;
#line 187 "csptri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 188 "csptri.f"
	    i__2 = kp;
#line 188 "csptri.f"
	    if (ipiv[*info] > 0 && (ap[i__2].r == 0. && ap[i__2].i == 0.)) {
#line 188 "csptri.f"
		return 0;
#line 188 "csptri.f"
	    }
#line 190 "csptri.f"
	    kp = kp + *n - *info + 1;
#line 191 "csptri.f"
/* L20: */
#line 191 "csptri.f"
	}
#line 192 "csptri.f"
    }
#line 193 "csptri.f"
    *info = 0;

#line 195 "csptri.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 202 "csptri.f"
	k = 1;
#line 203 "csptri.f"
	kc = 1;
#line 204 "csptri.f"
L30:

/*        If K > N, exit from loop. */

#line 208 "csptri.f"
	if (k > *n) {
#line 208 "csptri.f"
	    goto L50;
#line 208 "csptri.f"
	}

#line 211 "csptri.f"
	kcnext = kc + k;
#line 212 "csptri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 218 "csptri.f"
	    i__1 = kc + k - 1;
#line 218 "csptri.f"
	    z_div(&z__1, &c_b1, &ap[kc + k - 1]);
#line 218 "csptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;

/*           Compute column K of the inverse. */

#line 222 "csptri.f"
	    if (k > 1) {
#line 223 "csptri.f"
		i__1 = k - 1;
#line 223 "csptri.f"
		ccopy_(&i__1, &ap[kc], &c__1, &work[1], &c__1);
#line 224 "csptri.f"
		i__1 = k - 1;
#line 224 "csptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 224 "csptri.f"
		cspmv_(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &
			ap[kc], &c__1, (ftnlen)1);
#line 226 "csptri.f"
		i__1 = kc + k - 1;
#line 226 "csptri.f"
		i__2 = kc + k - 1;
#line 226 "csptri.f"
		i__3 = k - 1;
#line 226 "csptri.f"
		cdotu_(&z__2, &i__3, &work[1], &c__1, &ap[kc], &c__1);
#line 226 "csptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 226 "csptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 228 "csptri.f"
	    }
#line 229 "csptri.f"
	    kstep = 1;
#line 230 "csptri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 236 "csptri.f"
	    i__1 = kcnext + k - 1;
#line 236 "csptri.f"
	    t.r = ap[i__1].r, t.i = ap[i__1].i;
#line 237 "csptri.f"
	    z_div(&z__1, &ap[kc + k - 1], &t);
#line 237 "csptri.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 238 "csptri.f"
	    z_div(&z__1, &ap[kcnext + k], &t);
#line 238 "csptri.f"
	    akp1.r = z__1.r, akp1.i = z__1.i;
#line 239 "csptri.f"
	    z_div(&z__1, &ap[kcnext + k - 1], &t);
#line 239 "csptri.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 240 "csptri.f"
	    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + 
		    ak.i * akp1.r;
#line 240 "csptri.f"
	    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 240 "csptri.f"
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
#line 240 "csptri.f"
	    d__.r = z__1.r, d__.i = z__1.i;
#line 241 "csptri.f"
	    i__1 = kc + k - 1;
#line 241 "csptri.f"
	    z_div(&z__1, &akp1, &d__);
#line 241 "csptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 242 "csptri.f"
	    i__1 = kcnext + k;
#line 242 "csptri.f"
	    z_div(&z__1, &ak, &d__);
#line 242 "csptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 243 "csptri.f"
	    i__1 = kcnext + k - 1;
#line 243 "csptri.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 243 "csptri.f"
	    z_div(&z__1, &z__2, &d__);
#line 243 "csptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;

/*           Compute columns K and K+1 of the inverse. */

#line 247 "csptri.f"
	    if (k > 1) {
#line 248 "csptri.f"
		i__1 = k - 1;
#line 248 "csptri.f"
		ccopy_(&i__1, &ap[kc], &c__1, &work[1], &c__1);
#line 249 "csptri.f"
		i__1 = k - 1;
#line 249 "csptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 249 "csptri.f"
		cspmv_(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &
			ap[kc], &c__1, (ftnlen)1);
#line 251 "csptri.f"
		i__1 = kc + k - 1;
#line 251 "csptri.f"
		i__2 = kc + k - 1;
#line 251 "csptri.f"
		i__3 = k - 1;
#line 251 "csptri.f"
		cdotu_(&z__2, &i__3, &work[1], &c__1, &ap[kc], &c__1);
#line 251 "csptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 251 "csptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 253 "csptri.f"
		i__1 = kcnext + k - 1;
#line 253 "csptri.f"
		i__2 = kcnext + k - 1;
#line 253 "csptri.f"
		i__3 = k - 1;
#line 253 "csptri.f"
		cdotu_(&z__2, &i__3, &ap[kc], &c__1, &ap[kcnext], &c__1);
#line 253 "csptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 253 "csptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 256 "csptri.f"
		i__1 = k - 1;
#line 256 "csptri.f"
		ccopy_(&i__1, &ap[kcnext], &c__1, &work[1], &c__1);
#line 257 "csptri.f"
		i__1 = k - 1;
#line 257 "csptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 257 "csptri.f"
		cspmv_(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &
			ap[kcnext], &c__1, (ftnlen)1);
#line 259 "csptri.f"
		i__1 = kcnext + k;
#line 259 "csptri.f"
		i__2 = kcnext + k;
#line 259 "csptri.f"
		i__3 = k - 1;
#line 259 "csptri.f"
		cdotu_(&z__2, &i__3, &work[1], &c__1, &ap[kcnext], &c__1);
#line 259 "csptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 259 "csptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 261 "csptri.f"
	    }
#line 262 "csptri.f"
	    kstep = 2;
#line 263 "csptri.f"
	    kcnext = kcnext + k + 1;
#line 264 "csptri.f"
	}

#line 266 "csptri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 267 "csptri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 272 "csptri.f"
	    kpc = (kp - 1) * kp / 2 + 1;
#line 273 "csptri.f"
	    i__1 = kp - 1;
#line 273 "csptri.f"
	    cswap_(&i__1, &ap[kc], &c__1, &ap[kpc], &c__1);
#line 274 "csptri.f"
	    kx = kpc + kp - 1;
#line 275 "csptri.f"
	    i__1 = k - 1;
#line 275 "csptri.f"
	    for (j = kp + 1; j <= i__1; ++j) {
#line 276 "csptri.f"
		kx = kx + j - 1;
#line 277 "csptri.f"
		i__2 = kc + j - 1;
#line 277 "csptri.f"
		temp.r = ap[i__2].r, temp.i = ap[i__2].i;
#line 278 "csptri.f"
		i__2 = kc + j - 1;
#line 278 "csptri.f"
		i__3 = kx;
#line 278 "csptri.f"
		ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
#line 279 "csptri.f"
		i__2 = kx;
#line 279 "csptri.f"
		ap[i__2].r = temp.r, ap[i__2].i = temp.i;
#line 280 "csptri.f"
/* L40: */
#line 280 "csptri.f"
	    }
#line 281 "csptri.f"
	    i__1 = kc + k - 1;
#line 281 "csptri.f"
	    temp.r = ap[i__1].r, temp.i = ap[i__1].i;
#line 282 "csptri.f"
	    i__1 = kc + k - 1;
#line 282 "csptri.f"
	    i__2 = kpc + kp - 1;
#line 282 "csptri.f"
	    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 283 "csptri.f"
	    i__1 = kpc + kp - 1;
#line 283 "csptri.f"
	    ap[i__1].r = temp.r, ap[i__1].i = temp.i;
#line 284 "csptri.f"
	    if (kstep == 2) {
#line 285 "csptri.f"
		i__1 = kc + k + k - 1;
#line 285 "csptri.f"
		temp.r = ap[i__1].r, temp.i = ap[i__1].i;
#line 286 "csptri.f"
		i__1 = kc + k + k - 1;
#line 286 "csptri.f"
		i__2 = kc + k + kp - 1;
#line 286 "csptri.f"
		ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 287 "csptri.f"
		i__1 = kc + k + kp - 1;
#line 287 "csptri.f"
		ap[i__1].r = temp.r, ap[i__1].i = temp.i;
#line 288 "csptri.f"
	    }
#line 289 "csptri.f"
	}

#line 291 "csptri.f"
	k += kstep;
#line 292 "csptri.f"
	kc = kcnext;
#line 293 "csptri.f"
	goto L30;
#line 294 "csptri.f"
L50:

#line 296 "csptri.f"
	;
#line 296 "csptri.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 303 "csptri.f"
	npp = *n * (*n + 1) / 2;
#line 304 "csptri.f"
	k = *n;
#line 305 "csptri.f"
	kc = npp;
#line 306 "csptri.f"
L60:

/*        If K < 1, exit from loop. */

#line 310 "csptri.f"
	if (k < 1) {
#line 310 "csptri.f"
	    goto L80;
#line 310 "csptri.f"
	}

#line 313 "csptri.f"
	kcnext = kc - (*n - k + 2);
#line 314 "csptri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 320 "csptri.f"
	    i__1 = kc;
#line 320 "csptri.f"
	    z_div(&z__1, &c_b1, &ap[kc]);
#line 320 "csptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;

/*           Compute column K of the inverse. */

#line 324 "csptri.f"
	    if (k < *n) {
#line 325 "csptri.f"
		i__1 = *n - k;
#line 325 "csptri.f"
		ccopy_(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
#line 326 "csptri.f"
		i__1 = *n - k;
#line 326 "csptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 326 "csptri.f"
		cspmv_(uplo, &i__1, &z__1, &ap[kc + *n - k + 1], &work[1], &
			c__1, &c_b2, &ap[kc + 1], &c__1, (ftnlen)1);
#line 328 "csptri.f"
		i__1 = kc;
#line 328 "csptri.f"
		i__2 = kc;
#line 328 "csptri.f"
		i__3 = *n - k;
#line 328 "csptri.f"
		cdotu_(&z__2, &i__3, &work[1], &c__1, &ap[kc + 1], &c__1);
#line 328 "csptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 328 "csptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 330 "csptri.f"
	    }
#line 331 "csptri.f"
	    kstep = 1;
#line 332 "csptri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 338 "csptri.f"
	    i__1 = kcnext + 1;
#line 338 "csptri.f"
	    t.r = ap[i__1].r, t.i = ap[i__1].i;
#line 339 "csptri.f"
	    z_div(&z__1, &ap[kcnext], &t);
#line 339 "csptri.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 340 "csptri.f"
	    z_div(&z__1, &ap[kc], &t);
#line 340 "csptri.f"
	    akp1.r = z__1.r, akp1.i = z__1.i;
#line 341 "csptri.f"
	    z_div(&z__1, &ap[kcnext + 1], &t);
#line 341 "csptri.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 342 "csptri.f"
	    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + 
		    ak.i * akp1.r;
#line 342 "csptri.f"
	    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
#line 342 "csptri.f"
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
#line 342 "csptri.f"
	    d__.r = z__1.r, d__.i = z__1.i;
#line 343 "csptri.f"
	    i__1 = kcnext;
#line 343 "csptri.f"
	    z_div(&z__1, &akp1, &d__);
#line 343 "csptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 344 "csptri.f"
	    i__1 = kc;
#line 344 "csptri.f"
	    z_div(&z__1, &ak, &d__);
#line 344 "csptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 345 "csptri.f"
	    i__1 = kcnext + 1;
#line 345 "csptri.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 345 "csptri.f"
	    z_div(&z__1, &z__2, &d__);
#line 345 "csptri.f"
	    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;

/*           Compute columns K-1 and K of the inverse. */

#line 349 "csptri.f"
	    if (k < *n) {
#line 350 "csptri.f"
		i__1 = *n - k;
#line 350 "csptri.f"
		ccopy_(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
#line 351 "csptri.f"
		i__1 = *n - k;
#line 351 "csptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 351 "csptri.f"
		cspmv_(uplo, &i__1, &z__1, &ap[kc + (*n - k + 1)], &work[1], &
			c__1, &c_b2, &ap[kc + 1], &c__1, (ftnlen)1);
#line 353 "csptri.f"
		i__1 = kc;
#line 353 "csptri.f"
		i__2 = kc;
#line 353 "csptri.f"
		i__3 = *n - k;
#line 353 "csptri.f"
		cdotu_(&z__2, &i__3, &work[1], &c__1, &ap[kc + 1], &c__1);
#line 353 "csptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 353 "csptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 355 "csptri.f"
		i__1 = kcnext + 1;
#line 355 "csptri.f"
		i__2 = kcnext + 1;
#line 355 "csptri.f"
		i__3 = *n - k;
#line 355 "csptri.f"
		cdotu_(&z__2, &i__3, &ap[kc + 1], &c__1, &ap[kcnext + 2], &
			c__1);
#line 355 "csptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 355 "csptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 358 "csptri.f"
		i__1 = *n - k;
#line 358 "csptri.f"
		ccopy_(&i__1, &ap[kcnext + 2], &c__1, &work[1], &c__1);
#line 359 "csptri.f"
		i__1 = *n - k;
#line 359 "csptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 359 "csptri.f"
		cspmv_(uplo, &i__1, &z__1, &ap[kc + (*n - k + 1)], &work[1], &
			c__1, &c_b2, &ap[kcnext + 2], &c__1, (ftnlen)1);
#line 361 "csptri.f"
		i__1 = kcnext;
#line 361 "csptri.f"
		i__2 = kcnext;
#line 361 "csptri.f"
		i__3 = *n - k;
#line 361 "csptri.f"
		cdotu_(&z__2, &i__3, &work[1], &c__1, &ap[kcnext + 2], &c__1);
#line 361 "csptri.f"
		z__1.r = ap[i__2].r - z__2.r, z__1.i = ap[i__2].i - z__2.i;
#line 361 "csptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 363 "csptri.f"
	    }
#line 364 "csptri.f"
	    kstep = 2;
#line 365 "csptri.f"
	    kcnext -= *n - k + 3;
#line 366 "csptri.f"
	}

#line 368 "csptri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 369 "csptri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 374 "csptri.f"
	    kpc = npp - (*n - kp + 1) * (*n - kp + 2) / 2 + 1;
#line 375 "csptri.f"
	    if (kp < *n) {
#line 375 "csptri.f"
		i__1 = *n - kp;
#line 375 "csptri.f"
		cswap_(&i__1, &ap[kc + kp - k + 1], &c__1, &ap[kpc + 1], &
			c__1);
#line 375 "csptri.f"
	    }
#line 377 "csptri.f"
	    kx = kc + kp - k;
#line 378 "csptri.f"
	    i__1 = kp - 1;
#line 378 "csptri.f"
	    for (j = k + 1; j <= i__1; ++j) {
#line 379 "csptri.f"
		kx = kx + *n - j + 1;
#line 380 "csptri.f"
		i__2 = kc + j - k;
#line 380 "csptri.f"
		temp.r = ap[i__2].r, temp.i = ap[i__2].i;
#line 381 "csptri.f"
		i__2 = kc + j - k;
#line 381 "csptri.f"
		i__3 = kx;
#line 381 "csptri.f"
		ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
#line 382 "csptri.f"
		i__2 = kx;
#line 382 "csptri.f"
		ap[i__2].r = temp.r, ap[i__2].i = temp.i;
#line 383 "csptri.f"
/* L70: */
#line 383 "csptri.f"
	    }
#line 384 "csptri.f"
	    i__1 = kc;
#line 384 "csptri.f"
	    temp.r = ap[i__1].r, temp.i = ap[i__1].i;
#line 385 "csptri.f"
	    i__1 = kc;
#line 385 "csptri.f"
	    i__2 = kpc;
#line 385 "csptri.f"
	    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 386 "csptri.f"
	    i__1 = kpc;
#line 386 "csptri.f"
	    ap[i__1].r = temp.r, ap[i__1].i = temp.i;
#line 387 "csptri.f"
	    if (kstep == 2) {
#line 388 "csptri.f"
		i__1 = kc - *n + k - 1;
#line 388 "csptri.f"
		temp.r = ap[i__1].r, temp.i = ap[i__1].i;
#line 389 "csptri.f"
		i__1 = kc - *n + k - 1;
#line 389 "csptri.f"
		i__2 = kc - *n + kp - 1;
#line 389 "csptri.f"
		ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
#line 390 "csptri.f"
		i__1 = kc - *n + kp - 1;
#line 390 "csptri.f"
		ap[i__1].r = temp.r, ap[i__1].i = temp.i;
#line 391 "csptri.f"
	    }
#line 392 "csptri.f"
	}

#line 394 "csptri.f"
	k -= kstep;
#line 395 "csptri.f"
	kc = kcnext;
#line 396 "csptri.f"
	goto L60;
#line 397 "csptri.f"
L80:
#line 398 "csptri.f"
	;
#line 398 "csptri.f"
    }

#line 400 "csptri.f"
    return 0;

/*     End of CSPTRI */

} /* csptri_ */

