#line 1 "chetri_rook.f"
/* chetri_rook.f -- translated by f2c (version 20100827).
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

#line 1 "chetri_rook.f"
/* Table of constant values */

static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b CHETRI_ROOK computes the inverse of HE matrix using the factorization obtained with the bounded
 Bunch-Kaufman ("rook") diagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRI_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO ) */

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
/* > CHETRI_ROOK computes the inverse of a complex Hermitian indefinite matrix */
/* > A using the factorization A = U*D*U**H or A = L*D*L**H computed by */
/* > CHETRF_ROOK. */
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
/* >          used to obtain the factor U or L as computed by CHETRF_ROOK. */
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
/* >          as determined by CHETRF_ROOK. */
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

/* > \date November 2013 */

/* > \ingroup complexHEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2013,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int chetri_rook__(char *uplo, integer *n, doublecomplex *a, 
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


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

#line 174 "chetri_rook.f"
    /* Parameter adjustments */
#line 174 "chetri_rook.f"
    a_dim1 = *lda;
#line 174 "chetri_rook.f"
    a_offset = 1 + a_dim1;
#line 174 "chetri_rook.f"
    a -= a_offset;
#line 174 "chetri_rook.f"
    --ipiv;
#line 174 "chetri_rook.f"
    --work;
#line 174 "chetri_rook.f"

#line 174 "chetri_rook.f"
    /* Function Body */
#line 174 "chetri_rook.f"
    *info = 0;
#line 175 "chetri_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 176 "chetri_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 177 "chetri_rook.f"
	*info = -1;
#line 178 "chetri_rook.f"
    } else if (*n < 0) {
#line 179 "chetri_rook.f"
	*info = -2;
#line 180 "chetri_rook.f"
    } else if (*lda < max(1,*n)) {
#line 181 "chetri_rook.f"
	*info = -4;
#line 182 "chetri_rook.f"
    }
#line 183 "chetri_rook.f"
    if (*info != 0) {
#line 184 "chetri_rook.f"
	i__1 = -(*info);
#line 184 "chetri_rook.f"
	xerbla_("CHETRI_ROOK", &i__1, (ftnlen)11);
#line 185 "chetri_rook.f"
	return 0;
#line 186 "chetri_rook.f"
    }

/*     Quick return if possible */

#line 190 "chetri_rook.f"
    if (*n == 0) {
#line 190 "chetri_rook.f"
	return 0;
#line 190 "chetri_rook.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 195 "chetri_rook.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 199 "chetri_rook.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 200 "chetri_rook.f"
	    i__1 = *info + *info * a_dim1;
#line 200 "chetri_rook.f"
	    if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
#line 200 "chetri_rook.f"
		return 0;
#line 200 "chetri_rook.f"
	    }
#line 202 "chetri_rook.f"
/* L10: */
#line 202 "chetri_rook.f"
	}
#line 203 "chetri_rook.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 207 "chetri_rook.f"
	i__1 = *n;
#line 207 "chetri_rook.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 208 "chetri_rook.f"
	    i__2 = *info + *info * a_dim1;
#line 208 "chetri_rook.f"
	    if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
#line 208 "chetri_rook.f"
		return 0;
#line 208 "chetri_rook.f"
	    }
#line 210 "chetri_rook.f"
/* L20: */
#line 210 "chetri_rook.f"
	}
#line 211 "chetri_rook.f"
    }
#line 212 "chetri_rook.f"
    *info = 0;

#line 214 "chetri_rook.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**H. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 221 "chetri_rook.f"
	k = 1;
#line 222 "chetri_rook.f"
L30:

/*        If K > N, exit from loop. */

#line 226 "chetri_rook.f"
	if (k > *n) {
#line 226 "chetri_rook.f"
	    goto L70;
#line 226 "chetri_rook.f"
	}

#line 229 "chetri_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 235 "chetri_rook.f"
	    i__1 = k + k * a_dim1;
#line 235 "chetri_rook.f"
	    i__2 = k + k * a_dim1;
#line 235 "chetri_rook.f"
	    d__1 = 1. / a[i__2].r;
#line 235 "chetri_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;

/*           Compute column K of the inverse. */

#line 239 "chetri_rook.f"
	    if (k > 1) {
#line 240 "chetri_rook.f"
		i__1 = k - 1;
#line 240 "chetri_rook.f"
		ccopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 241 "chetri_rook.f"
		i__1 = k - 1;
#line 241 "chetri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 241 "chetri_rook.f"
		chemv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 243 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 243 "chetri_rook.f"
		i__2 = k + k * a_dim1;
#line 243 "chetri_rook.f"
		i__3 = k - 1;
#line 243 "chetri_rook.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 243 "chetri_rook.f"
		d__1 = z__2.r;
#line 243 "chetri_rook.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 243 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 245 "chetri_rook.f"
	    }
#line 246 "chetri_rook.f"
	    kstep = 1;
#line 247 "chetri_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 253 "chetri_rook.f"
	    t = z_abs(&a[k + (k + 1) * a_dim1]);
#line 254 "chetri_rook.f"
	    i__1 = k + k * a_dim1;
#line 254 "chetri_rook.f"
	    ak = a[i__1].r / t;
#line 255 "chetri_rook.f"
	    i__1 = k + 1 + (k + 1) * a_dim1;
#line 255 "chetri_rook.f"
	    akp1 = a[i__1].r / t;
#line 256 "chetri_rook.f"
	    i__1 = k + (k + 1) * a_dim1;
#line 256 "chetri_rook.f"
	    z__1.r = a[i__1].r / t, z__1.i = a[i__1].i / t;
#line 256 "chetri_rook.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 257 "chetri_rook.f"
	    d__ = t * (ak * akp1 - 1.);
#line 258 "chetri_rook.f"
	    i__1 = k + k * a_dim1;
#line 258 "chetri_rook.f"
	    d__1 = akp1 / d__;
#line 258 "chetri_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 259 "chetri_rook.f"
	    i__1 = k + 1 + (k + 1) * a_dim1;
#line 259 "chetri_rook.f"
	    d__1 = ak / d__;
#line 259 "chetri_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 260 "chetri_rook.f"
	    i__1 = k + (k + 1) * a_dim1;
#line 260 "chetri_rook.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 260 "chetri_rook.f"
	    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
#line 260 "chetri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute columns K and K+1 of the inverse. */

#line 264 "chetri_rook.f"
	    if (k > 1) {
#line 265 "chetri_rook.f"
		i__1 = k - 1;
#line 265 "chetri_rook.f"
		ccopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 266 "chetri_rook.f"
		i__1 = k - 1;
#line 266 "chetri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 266 "chetri_rook.f"
		chemv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 268 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 268 "chetri_rook.f"
		i__2 = k + k * a_dim1;
#line 268 "chetri_rook.f"
		i__3 = k - 1;
#line 268 "chetri_rook.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
#line 268 "chetri_rook.f"
		d__1 = z__2.r;
#line 268 "chetri_rook.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 268 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 270 "chetri_rook.f"
		i__1 = k + (k + 1) * a_dim1;
#line 270 "chetri_rook.f"
		i__2 = k + (k + 1) * a_dim1;
#line 270 "chetri_rook.f"
		i__3 = k - 1;
#line 270 "chetri_rook.f"
		cdotc_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * 
			a_dim1 + 1], &c__1);
#line 270 "chetri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 270 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 272 "chetri_rook.f"
		i__1 = k - 1;
#line 272 "chetri_rook.f"
		ccopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &
			c__1);
#line 273 "chetri_rook.f"
		i__1 = k - 1;
#line 273 "chetri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 273 "chetri_rook.f"
		chemv_(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1,
			 &c_b2, &a[(k + 1) * a_dim1 + 1], &c__1, (ftnlen)1);
#line 275 "chetri_rook.f"
		i__1 = k + 1 + (k + 1) * a_dim1;
#line 275 "chetri_rook.f"
		i__2 = k + 1 + (k + 1) * a_dim1;
#line 275 "chetri_rook.f"
		i__3 = k - 1;
#line 275 "chetri_rook.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1]
			, &c__1);
#line 275 "chetri_rook.f"
		d__1 = z__2.r;
#line 275 "chetri_rook.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 275 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 278 "chetri_rook.f"
	    }
#line 279 "chetri_rook.f"
	    kstep = 2;
#line 280 "chetri_rook.f"
	}

#line 282 "chetri_rook.f"
	if (kstep == 1) {

/*           Interchange rows and columns K and IPIV(K) in the leading */
/*           submatrix A(1:k,1:k) */

#line 287 "chetri_rook.f"
	    kp = ipiv[k];
#line 288 "chetri_rook.f"
	    if (kp != k) {

#line 290 "chetri_rook.f"
		if (kp > 1) {
#line 290 "chetri_rook.f"
		    i__1 = kp - 1;
#line 290 "chetri_rook.f"
		    cswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 290 "chetri_rook.f"
		}

#line 293 "chetri_rook.f"
		i__1 = k - 1;
#line 293 "chetri_rook.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 294 "chetri_rook.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 294 "chetri_rook.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 295 "chetri_rook.f"
		    i__2 = j + k * a_dim1;
#line 295 "chetri_rook.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 295 "chetri_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 296 "chetri_rook.f"
		    i__2 = kp + j * a_dim1;
#line 296 "chetri_rook.f"
		    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 297 "chetri_rook.f"
/* L40: */
#line 297 "chetri_rook.f"
		}

#line 299 "chetri_rook.f"
		i__1 = kp + k * a_dim1;
#line 299 "chetri_rook.f"
		d_cnjg(&z__1, &a[kp + k * a_dim1]);
#line 299 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 301 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 301 "chetri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 302 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 302 "chetri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 302 "chetri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 303 "chetri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 303 "chetri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 304 "chetri_rook.f"
	    }
#line 305 "chetri_rook.f"
	} else {

/*           Interchange rows and columns K and K+1 with -IPIV(K) and */
/*           -IPIV(K+1) in the leading submatrix A(k+1:n,k+1:n) */

/*           (1) Interchange rows and columns K and -IPIV(K) */

#line 312 "chetri_rook.f"
	    kp = -ipiv[k];
#line 313 "chetri_rook.f"
	    if (kp != k) {

#line 315 "chetri_rook.f"
		if (kp > 1) {
#line 315 "chetri_rook.f"
		    i__1 = kp - 1;
#line 315 "chetri_rook.f"
		    cswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 315 "chetri_rook.f"
		}

#line 318 "chetri_rook.f"
		i__1 = k - 1;
#line 318 "chetri_rook.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 319 "chetri_rook.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 319 "chetri_rook.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 320 "chetri_rook.f"
		    i__2 = j + k * a_dim1;
#line 320 "chetri_rook.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 320 "chetri_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 321 "chetri_rook.f"
		    i__2 = kp + j * a_dim1;
#line 321 "chetri_rook.f"
		    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 322 "chetri_rook.f"
/* L50: */
#line 322 "chetri_rook.f"
		}

#line 324 "chetri_rook.f"
		i__1 = kp + k * a_dim1;
#line 324 "chetri_rook.f"
		d_cnjg(&z__1, &a[kp + k * a_dim1]);
#line 324 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 326 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 326 "chetri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 327 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 327 "chetri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 327 "chetri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 328 "chetri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 328 "chetri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;

#line 330 "chetri_rook.f"
		i__1 = k + (k + 1) * a_dim1;
#line 330 "chetri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 331 "chetri_rook.f"
		i__1 = k + (k + 1) * a_dim1;
#line 331 "chetri_rook.f"
		i__2 = kp + (k + 1) * a_dim1;
#line 331 "chetri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 332 "chetri_rook.f"
		i__1 = kp + (k + 1) * a_dim1;
#line 332 "chetri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 333 "chetri_rook.f"
	    }

/*           (2) Interchange rows and columns K+1 and -IPIV(K+1) */

#line 337 "chetri_rook.f"
	    ++k;
#line 338 "chetri_rook.f"
	    kp = -ipiv[k];
#line 339 "chetri_rook.f"
	    if (kp != k) {

#line 341 "chetri_rook.f"
		if (kp > 1) {
#line 341 "chetri_rook.f"
		    i__1 = kp - 1;
#line 341 "chetri_rook.f"
		    cswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 341 "chetri_rook.f"
		}

#line 344 "chetri_rook.f"
		i__1 = k - 1;
#line 344 "chetri_rook.f"
		for (j = kp + 1; j <= i__1; ++j) {
#line 345 "chetri_rook.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 345 "chetri_rook.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 346 "chetri_rook.f"
		    i__2 = j + k * a_dim1;
#line 346 "chetri_rook.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 346 "chetri_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 347 "chetri_rook.f"
		    i__2 = kp + j * a_dim1;
#line 347 "chetri_rook.f"
		    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 348 "chetri_rook.f"
/* L60: */
#line 348 "chetri_rook.f"
		}

#line 350 "chetri_rook.f"
		i__1 = kp + k * a_dim1;
#line 350 "chetri_rook.f"
		d_cnjg(&z__1, &a[kp + k * a_dim1]);
#line 350 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 352 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 352 "chetri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 353 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 353 "chetri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 353 "chetri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 354 "chetri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 354 "chetri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 355 "chetri_rook.f"
	    }
#line 356 "chetri_rook.f"
	}

#line 358 "chetri_rook.f"
	++k;
#line 359 "chetri_rook.f"
	goto L30;
#line 360 "chetri_rook.f"
L70:

#line 362 "chetri_rook.f"
	;
#line 362 "chetri_rook.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**H. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 369 "chetri_rook.f"
	k = *n;
#line 370 "chetri_rook.f"
L80:

/*        If K < 1, exit from loop. */

#line 374 "chetri_rook.f"
	if (k < 1) {
#line 374 "chetri_rook.f"
	    goto L120;
#line 374 "chetri_rook.f"
	}

#line 377 "chetri_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 383 "chetri_rook.f"
	    i__1 = k + k * a_dim1;
#line 383 "chetri_rook.f"
	    i__2 = k + k * a_dim1;
#line 383 "chetri_rook.f"
	    d__1 = 1. / a[i__2].r;
#line 383 "chetri_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;

/*           Compute column K of the inverse. */

#line 387 "chetri_rook.f"
	    if (k < *n) {
#line 388 "chetri_rook.f"
		i__1 = *n - k;
#line 388 "chetri_rook.f"
		ccopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 389 "chetri_rook.f"
		i__1 = *n - k;
#line 389 "chetri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 389 "chetri_rook.f"
		chemv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1,
			 (ftnlen)1);
#line 391 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 391 "chetri_rook.f"
		i__2 = k + k * a_dim1;
#line 391 "chetri_rook.f"
		i__3 = *n - k;
#line 391 "chetri_rook.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], 
			&c__1);
#line 391 "chetri_rook.f"
		d__1 = z__2.r;
#line 391 "chetri_rook.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 391 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 393 "chetri_rook.f"
	    }
#line 394 "chetri_rook.f"
	    kstep = 1;
#line 395 "chetri_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 401 "chetri_rook.f"
	    t = z_abs(&a[k + (k - 1) * a_dim1]);
#line 402 "chetri_rook.f"
	    i__1 = k - 1 + (k - 1) * a_dim1;
#line 402 "chetri_rook.f"
	    ak = a[i__1].r / t;
#line 403 "chetri_rook.f"
	    i__1 = k + k * a_dim1;
#line 403 "chetri_rook.f"
	    akp1 = a[i__1].r / t;
#line 404 "chetri_rook.f"
	    i__1 = k + (k - 1) * a_dim1;
#line 404 "chetri_rook.f"
	    z__1.r = a[i__1].r / t, z__1.i = a[i__1].i / t;
#line 404 "chetri_rook.f"
	    akkp1.r = z__1.r, akkp1.i = z__1.i;
#line 405 "chetri_rook.f"
	    d__ = t * (ak * akp1 - 1.);
#line 406 "chetri_rook.f"
	    i__1 = k - 1 + (k - 1) * a_dim1;
#line 406 "chetri_rook.f"
	    d__1 = akp1 / d__;
#line 406 "chetri_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 407 "chetri_rook.f"
	    i__1 = k + k * a_dim1;
#line 407 "chetri_rook.f"
	    d__1 = ak / d__;
#line 407 "chetri_rook.f"
	    a[i__1].r = d__1, a[i__1].i = 0.;
#line 408 "chetri_rook.f"
	    i__1 = k + (k - 1) * a_dim1;
#line 408 "chetri_rook.f"
	    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
#line 408 "chetri_rook.f"
	    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
#line 408 "chetri_rook.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*           Compute columns K-1 and K of the inverse. */

#line 412 "chetri_rook.f"
	    if (k < *n) {
#line 413 "chetri_rook.f"
		i__1 = *n - k;
#line 413 "chetri_rook.f"
		ccopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 414 "chetri_rook.f"
		i__1 = *n - k;
#line 414 "chetri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 414 "chetri_rook.f"
		chemv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1,
			 (ftnlen)1);
#line 416 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 416 "chetri_rook.f"
		i__2 = k + k * a_dim1;
#line 416 "chetri_rook.f"
		i__3 = *n - k;
#line 416 "chetri_rook.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], 
			&c__1);
#line 416 "chetri_rook.f"
		d__1 = z__2.r;
#line 416 "chetri_rook.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 416 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 418 "chetri_rook.f"
		i__1 = k + (k - 1) * a_dim1;
#line 418 "chetri_rook.f"
		i__2 = k + (k - 1) * a_dim1;
#line 418 "chetri_rook.f"
		i__3 = *n - k;
#line 418 "chetri_rook.f"
		cdotc_(&z__2, &i__3, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 
			+ (k - 1) * a_dim1], &c__1);
#line 418 "chetri_rook.f"
		z__1.r = a[i__2].r - z__2.r, z__1.i = a[i__2].i - z__2.i;
#line 418 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 421 "chetri_rook.f"
		i__1 = *n - k;
#line 421 "chetri_rook.f"
		ccopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &
			c__1);
#line 422 "chetri_rook.f"
		i__1 = *n - k;
#line 422 "chetri_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 422 "chetri_rook.f"
		chemv_(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, 
			&work[1], &c__1, &c_b2, &a[k + 1 + (k - 1) * a_dim1], 
			&c__1, (ftnlen)1);
#line 424 "chetri_rook.f"
		i__1 = k - 1 + (k - 1) * a_dim1;
#line 424 "chetri_rook.f"
		i__2 = k - 1 + (k - 1) * a_dim1;
#line 424 "chetri_rook.f"
		i__3 = *n - k;
#line 424 "chetri_rook.f"
		cdotc_(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + (k - 1) * 
			a_dim1], &c__1);
#line 424 "chetri_rook.f"
		d__1 = z__2.r;
#line 424 "chetri_rook.f"
		z__1.r = a[i__2].r - d__1, z__1.i = a[i__2].i;
#line 424 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 427 "chetri_rook.f"
	    }
#line 428 "chetri_rook.f"
	    kstep = 2;
#line 429 "chetri_rook.f"
	}

#line 431 "chetri_rook.f"
	if (kstep == 1) {

/*           Interchange rows and columns K and IPIV(K) in the trailing */
/*           submatrix A(k:n,k:n) */

#line 436 "chetri_rook.f"
	    kp = ipiv[k];
#line 437 "chetri_rook.f"
	    if (kp != k) {

#line 439 "chetri_rook.f"
		if (kp < *n) {
#line 439 "chetri_rook.f"
		    i__1 = *n - kp;
#line 439 "chetri_rook.f"
		    cswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 439 "chetri_rook.f"
		}

#line 442 "chetri_rook.f"
		i__1 = kp - 1;
#line 442 "chetri_rook.f"
		for (j = k + 1; j <= i__1; ++j) {
#line 443 "chetri_rook.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 443 "chetri_rook.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 444 "chetri_rook.f"
		    i__2 = j + k * a_dim1;
#line 444 "chetri_rook.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 444 "chetri_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 445 "chetri_rook.f"
		    i__2 = kp + j * a_dim1;
#line 445 "chetri_rook.f"
		    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 446 "chetri_rook.f"
/* L90: */
#line 446 "chetri_rook.f"
		}

#line 448 "chetri_rook.f"
		i__1 = kp + k * a_dim1;
#line 448 "chetri_rook.f"
		d_cnjg(&z__1, &a[kp + k * a_dim1]);
#line 448 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 450 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 450 "chetri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 451 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 451 "chetri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 451 "chetri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 452 "chetri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 452 "chetri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 453 "chetri_rook.f"
	    }
#line 454 "chetri_rook.f"
	} else {

/*           Interchange rows and columns K and K-1 with -IPIV(K) and */
/*           -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n) */

/*           (1) Interchange rows and columns K and -IPIV(K) */

#line 461 "chetri_rook.f"
	    kp = -ipiv[k];
#line 462 "chetri_rook.f"
	    if (kp != k) {

#line 464 "chetri_rook.f"
		if (kp < *n) {
#line 464 "chetri_rook.f"
		    i__1 = *n - kp;
#line 464 "chetri_rook.f"
		    cswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 464 "chetri_rook.f"
		}

#line 467 "chetri_rook.f"
		i__1 = kp - 1;
#line 467 "chetri_rook.f"
		for (j = k + 1; j <= i__1; ++j) {
#line 468 "chetri_rook.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 468 "chetri_rook.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 469 "chetri_rook.f"
		    i__2 = j + k * a_dim1;
#line 469 "chetri_rook.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 469 "chetri_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 470 "chetri_rook.f"
		    i__2 = kp + j * a_dim1;
#line 470 "chetri_rook.f"
		    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 471 "chetri_rook.f"
/* L100: */
#line 471 "chetri_rook.f"
		}

#line 473 "chetri_rook.f"
		i__1 = kp + k * a_dim1;
#line 473 "chetri_rook.f"
		d_cnjg(&z__1, &a[kp + k * a_dim1]);
#line 473 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 475 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 475 "chetri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 476 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 476 "chetri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 476 "chetri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 477 "chetri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 477 "chetri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;

#line 479 "chetri_rook.f"
		i__1 = k + (k - 1) * a_dim1;
#line 479 "chetri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 480 "chetri_rook.f"
		i__1 = k + (k - 1) * a_dim1;
#line 480 "chetri_rook.f"
		i__2 = kp + (k - 1) * a_dim1;
#line 480 "chetri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 481 "chetri_rook.f"
		i__1 = kp + (k - 1) * a_dim1;
#line 481 "chetri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 482 "chetri_rook.f"
	    }

/*           (2) Interchange rows and columns K-1 and -IPIV(K-1) */

#line 486 "chetri_rook.f"
	    --k;
#line 487 "chetri_rook.f"
	    kp = -ipiv[k];
#line 488 "chetri_rook.f"
	    if (kp != k) {

#line 490 "chetri_rook.f"
		if (kp < *n) {
#line 490 "chetri_rook.f"
		    i__1 = *n - kp;
#line 490 "chetri_rook.f"
		    cswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 490 "chetri_rook.f"
		}

#line 493 "chetri_rook.f"
		i__1 = kp - 1;
#line 493 "chetri_rook.f"
		for (j = k + 1; j <= i__1; ++j) {
#line 494 "chetri_rook.f"
		    d_cnjg(&z__1, &a[j + k * a_dim1]);
#line 494 "chetri_rook.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 495 "chetri_rook.f"
		    i__2 = j + k * a_dim1;
#line 495 "chetri_rook.f"
		    d_cnjg(&z__1, &a[kp + j * a_dim1]);
#line 495 "chetri_rook.f"
		    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 496 "chetri_rook.f"
		    i__2 = kp + j * a_dim1;
#line 496 "chetri_rook.f"
		    a[i__2].r = temp.r, a[i__2].i = temp.i;
#line 497 "chetri_rook.f"
/* L110: */
#line 497 "chetri_rook.f"
		}

#line 499 "chetri_rook.f"
		i__1 = kp + k * a_dim1;
#line 499 "chetri_rook.f"
		d_cnjg(&z__1, &a[kp + k * a_dim1]);
#line 499 "chetri_rook.f"
		a[i__1].r = z__1.r, a[i__1].i = z__1.i;

#line 501 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 501 "chetri_rook.f"
		temp.r = a[i__1].r, temp.i = a[i__1].i;
#line 502 "chetri_rook.f"
		i__1 = k + k * a_dim1;
#line 502 "chetri_rook.f"
		i__2 = kp + kp * a_dim1;
#line 502 "chetri_rook.f"
		a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
#line 503 "chetri_rook.f"
		i__1 = kp + kp * a_dim1;
#line 503 "chetri_rook.f"
		a[i__1].r = temp.r, a[i__1].i = temp.i;
#line 504 "chetri_rook.f"
	    }
#line 505 "chetri_rook.f"
	}

#line 507 "chetri_rook.f"
	--k;
#line 508 "chetri_rook.f"
	goto L80;
#line 509 "chetri_rook.f"
L120:
#line 510 "chetri_rook.f"
	;
#line 510 "chetri_rook.f"
    }

#line 512 "chetri_rook.f"
    return 0;

/*     End of CHETRI_ROOK */

} /* chetri_rook__ */

