#line 1 "zpbstf.f"
/* zpbstf.f -- translated by f2c (version 20100827).
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

#line 1 "zpbstf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = -1.;

/* > \brief \b ZPBSTF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPBSTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpbstf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpbstf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpbstf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPBSTF( UPLO, N, KD, AB, LDAB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPBSTF computes a split Cholesky factorization of a complex */
/* > Hermitian positive definite band matrix A. */
/* > */
/* > This routine is designed to be used in conjunction with ZHBGST. */
/* > */
/* > The factorization has the form  A = S**H*S  where S is a band matrix */
/* > of the same bandwidth as A and the following structure: */
/* > */
/* >   S = ( U    ) */
/* >       ( M  L ) */
/* > */
/* > where U is upper triangular of order m = (n+kd)/2, and L is lower */
/* > triangular of order n-m. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first kd+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, if INFO = 0, the factor S from the split Cholesky */
/* >          factorization A = S**H*S. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, the factorization could not be completed, */
/* >               because the updated element a(i,i) was negative; the */
/* >               matrix A is not positive definite. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The band storage scheme is illustrated by the following example, when */
/* >  N = 7, KD = 2: */
/* > */
/* >  S = ( s11  s12  s13                     ) */
/* >      (      s22  s23  s24                ) */
/* >      (           s33  s34                ) */
/* >      (                s44                ) */
/* >      (           s53  s54  s55           ) */
/* >      (                s64  s65  s66      ) */
/* >      (                     s75  s76  s77 ) */
/* > */
/* >  If UPLO = 'U', the array AB holds: */
/* > */
/* >  on entry:                          on exit: */
/* > */
/* >   *    *   a13  a24  a35  a46  a57   *    *   s13  s24  s53**H s64**H s75**H */
/* >   *   a12  a23  a34  a45  a56  a67   *   s12  s23  s34  s54**H s65**H s76**H */
/* >  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55    s66    s77 */
/* > */
/* >  If UPLO = 'L', the array AB holds: */
/* > */
/* >  on entry:                          on exit: */
/* > */
/* >  a11  a22  a33  a44  a55  a66  a77  s11    s22    s33    s44  s55  s66  s77 */
/* >  a21  a32  a43  a54  a65  a76   *   s12**H s23**H s34**H s54  s65  s76   * */
/* >  a31  a42  a53  a64  a64   *    *   s13**H s24**H s53    s64  s75   *    * */
/* > */
/* >  Array elements marked * are not used by the routine; s12**H denotes */
/* >  conjg(s12); the diagonal elements of S are real. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zpbstf_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, m, km;
    static doublereal ajj;
    static integer kld;
    extern /* Subroutine */ int zher_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *), zlacgv_(
	    integer *, doublecomplex *, integer *);


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

#line 194 "zpbstf.f"
    /* Parameter adjustments */
#line 194 "zpbstf.f"
    ab_dim1 = *ldab;
#line 194 "zpbstf.f"
    ab_offset = 1 + ab_dim1;
#line 194 "zpbstf.f"
    ab -= ab_offset;
#line 194 "zpbstf.f"

#line 194 "zpbstf.f"
    /* Function Body */
#line 194 "zpbstf.f"
    *info = 0;
#line 195 "zpbstf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 196 "zpbstf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 197 "zpbstf.f"
	*info = -1;
#line 198 "zpbstf.f"
    } else if (*n < 0) {
#line 199 "zpbstf.f"
	*info = -2;
#line 200 "zpbstf.f"
    } else if (*kd < 0) {
#line 201 "zpbstf.f"
	*info = -3;
#line 202 "zpbstf.f"
    } else if (*ldab < *kd + 1) {
#line 203 "zpbstf.f"
	*info = -5;
#line 204 "zpbstf.f"
    }
#line 205 "zpbstf.f"
    if (*info != 0) {
#line 206 "zpbstf.f"
	i__1 = -(*info);
#line 206 "zpbstf.f"
	xerbla_("ZPBSTF", &i__1, (ftnlen)6);
#line 207 "zpbstf.f"
	return 0;
#line 208 "zpbstf.f"
    }

/*     Quick return if possible */

#line 212 "zpbstf.f"
    if (*n == 0) {
#line 212 "zpbstf.f"
	return 0;
#line 212 "zpbstf.f"
    }

/* Computing MAX */
#line 215 "zpbstf.f"
    i__1 = 1, i__2 = *ldab - 1;
#line 215 "zpbstf.f"
    kld = max(i__1,i__2);

/*     Set the splitting point m. */

#line 219 "zpbstf.f"
    m = (*n + *kd) / 2;

#line 221 "zpbstf.f"
    if (upper) {

/*        Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m). */

#line 225 "zpbstf.f"
	i__1 = m + 1;
#line 225 "zpbstf.f"
	for (j = *n; j >= i__1; --j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 229 "zpbstf.f"
	    i__2 = *kd + 1 + j * ab_dim1;
#line 229 "zpbstf.f"
	    ajj = ab[i__2].r;
#line 230 "zpbstf.f"
	    if (ajj <= 0.) {
#line 231 "zpbstf.f"
		i__2 = *kd + 1 + j * ab_dim1;
#line 231 "zpbstf.f"
		ab[i__2].r = ajj, ab[i__2].i = 0.;
#line 232 "zpbstf.f"
		goto L50;
#line 233 "zpbstf.f"
	    }
#line 234 "zpbstf.f"
	    ajj = sqrt(ajj);
#line 235 "zpbstf.f"
	    i__2 = *kd + 1 + j * ab_dim1;
#line 235 "zpbstf.f"
	    ab[i__2].r = ajj, ab[i__2].i = 0.;
/* Computing MIN */
#line 236 "zpbstf.f"
	    i__2 = j - 1;
#line 236 "zpbstf.f"
	    km = min(i__2,*kd);

/*           Compute elements j-km:j-1 of the j-th column and update the */
/*           the leading submatrix within the band. */

#line 241 "zpbstf.f"
	    d__1 = 1. / ajj;
#line 241 "zpbstf.f"
	    zdscal_(&km, &d__1, &ab[*kd + 1 - km + j * ab_dim1], &c__1);
#line 242 "zpbstf.f"
	    zher_("Upper", &km, &c_b9, &ab[*kd + 1 - km + j * ab_dim1], &c__1,
		     &ab[*kd + 1 + (j - km) * ab_dim1], &kld, (ftnlen)5);
#line 244 "zpbstf.f"
/* L10: */
#line 244 "zpbstf.f"
	}

/*        Factorize the updated submatrix A(1:m,1:m) as U**H*U. */

#line 248 "zpbstf.f"
	i__1 = m;
#line 248 "zpbstf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 252 "zpbstf.f"
	    i__2 = *kd + 1 + j * ab_dim1;
#line 252 "zpbstf.f"
	    ajj = ab[i__2].r;
#line 253 "zpbstf.f"
	    if (ajj <= 0.) {
#line 254 "zpbstf.f"
		i__2 = *kd + 1 + j * ab_dim1;
#line 254 "zpbstf.f"
		ab[i__2].r = ajj, ab[i__2].i = 0.;
#line 255 "zpbstf.f"
		goto L50;
#line 256 "zpbstf.f"
	    }
#line 257 "zpbstf.f"
	    ajj = sqrt(ajj);
#line 258 "zpbstf.f"
	    i__2 = *kd + 1 + j * ab_dim1;
#line 258 "zpbstf.f"
	    ab[i__2].r = ajj, ab[i__2].i = 0.;
/* Computing MIN */
#line 259 "zpbstf.f"
	    i__2 = *kd, i__3 = m - j;
#line 259 "zpbstf.f"
	    km = min(i__2,i__3);

/*           Compute elements j+1:j+km of the j-th row and update the */
/*           trailing submatrix within the band. */

#line 264 "zpbstf.f"
	    if (km > 0) {
#line 265 "zpbstf.f"
		d__1 = 1. / ajj;
#line 265 "zpbstf.f"
		zdscal_(&km, &d__1, &ab[*kd + (j + 1) * ab_dim1], &kld);
#line 266 "zpbstf.f"
		zlacgv_(&km, &ab[*kd + (j + 1) * ab_dim1], &kld);
#line 267 "zpbstf.f"
		zher_("Upper", &km, &c_b9, &ab[*kd + (j + 1) * ab_dim1], &kld,
			 &ab[*kd + 1 + (j + 1) * ab_dim1], &kld, (ftnlen)5);
#line 269 "zpbstf.f"
		zlacgv_(&km, &ab[*kd + (j + 1) * ab_dim1], &kld);
#line 270 "zpbstf.f"
	    }
#line 271 "zpbstf.f"
/* L20: */
#line 271 "zpbstf.f"
	}
#line 272 "zpbstf.f"
    } else {

/*        Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m). */

#line 276 "zpbstf.f"
	i__1 = m + 1;
#line 276 "zpbstf.f"
	for (j = *n; j >= i__1; --j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 280 "zpbstf.f"
	    i__2 = j * ab_dim1 + 1;
#line 280 "zpbstf.f"
	    ajj = ab[i__2].r;
#line 281 "zpbstf.f"
	    if (ajj <= 0.) {
#line 282 "zpbstf.f"
		i__2 = j * ab_dim1 + 1;
#line 282 "zpbstf.f"
		ab[i__2].r = ajj, ab[i__2].i = 0.;
#line 283 "zpbstf.f"
		goto L50;
#line 284 "zpbstf.f"
	    }
#line 285 "zpbstf.f"
	    ajj = sqrt(ajj);
#line 286 "zpbstf.f"
	    i__2 = j * ab_dim1 + 1;
#line 286 "zpbstf.f"
	    ab[i__2].r = ajj, ab[i__2].i = 0.;
/* Computing MIN */
#line 287 "zpbstf.f"
	    i__2 = j - 1;
#line 287 "zpbstf.f"
	    km = min(i__2,*kd);

/*           Compute elements j-km:j-1 of the j-th row and update the */
/*           trailing submatrix within the band. */

#line 292 "zpbstf.f"
	    d__1 = 1. / ajj;
#line 292 "zpbstf.f"
	    zdscal_(&km, &d__1, &ab[km + 1 + (j - km) * ab_dim1], &kld);
#line 293 "zpbstf.f"
	    zlacgv_(&km, &ab[km + 1 + (j - km) * ab_dim1], &kld);
#line 294 "zpbstf.f"
	    zher_("Lower", &km, &c_b9, &ab[km + 1 + (j - km) * ab_dim1], &kld,
		     &ab[(j - km) * ab_dim1 + 1], &kld, (ftnlen)5);
#line 296 "zpbstf.f"
	    zlacgv_(&km, &ab[km + 1 + (j - km) * ab_dim1], &kld);
#line 297 "zpbstf.f"
/* L30: */
#line 297 "zpbstf.f"
	}

/*        Factorize the updated submatrix A(1:m,1:m) as U**H*U. */

#line 301 "zpbstf.f"
	i__1 = m;
#line 301 "zpbstf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 305 "zpbstf.f"
	    i__2 = j * ab_dim1 + 1;
#line 305 "zpbstf.f"
	    ajj = ab[i__2].r;
#line 306 "zpbstf.f"
	    if (ajj <= 0.) {
#line 307 "zpbstf.f"
		i__2 = j * ab_dim1 + 1;
#line 307 "zpbstf.f"
		ab[i__2].r = ajj, ab[i__2].i = 0.;
#line 308 "zpbstf.f"
		goto L50;
#line 309 "zpbstf.f"
	    }
#line 310 "zpbstf.f"
	    ajj = sqrt(ajj);
#line 311 "zpbstf.f"
	    i__2 = j * ab_dim1 + 1;
#line 311 "zpbstf.f"
	    ab[i__2].r = ajj, ab[i__2].i = 0.;
/* Computing MIN */
#line 312 "zpbstf.f"
	    i__2 = *kd, i__3 = m - j;
#line 312 "zpbstf.f"
	    km = min(i__2,i__3);

/*           Compute elements j+1:j+km of the j-th column and update the */
/*           trailing submatrix within the band. */

#line 317 "zpbstf.f"
	    if (km > 0) {
#line 318 "zpbstf.f"
		d__1 = 1. / ajj;
#line 318 "zpbstf.f"
		zdscal_(&km, &d__1, &ab[j * ab_dim1 + 2], &c__1);
#line 319 "zpbstf.f"
		zher_("Lower", &km, &c_b9, &ab[j * ab_dim1 + 2], &c__1, &ab[(
			j + 1) * ab_dim1 + 1], &kld, (ftnlen)5);
#line 321 "zpbstf.f"
	    }
#line 322 "zpbstf.f"
/* L40: */
#line 322 "zpbstf.f"
	}
#line 323 "zpbstf.f"
    }
#line 324 "zpbstf.f"
    return 0;

#line 326 "zpbstf.f"
L50:
#line 327 "zpbstf.f"
    *info = j;
#line 328 "zpbstf.f"
    return 0;

/*     End of ZPBSTF */

} /* zpbstf_ */

