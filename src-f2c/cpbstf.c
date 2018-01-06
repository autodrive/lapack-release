#line 1 "cpbstf.f"
/* cpbstf.f -- translated by f2c (version 20100827).
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

#line 1 "cpbstf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = -1.;

/* > \brief \b CPBSTF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPBSTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbstf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbstf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbstf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPBSTF( UPLO, N, KD, AB, LDAB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPBSTF computes a split Cholesky factorization of a complex */
/* > Hermitian positive definite band matrix A. */
/* > */
/* > This routine is designed to be used in conjunction with CHBGST. */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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

/* > \ingroup complexOTHERcomputational */

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
/* Subroutine */ int cpbstf_(char *uplo, integer *n, integer *kd, 
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
    extern /* Subroutine */ int cher_(char *, integer *, doublereal *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , csscal_(integer *, doublereal *, doublecomplex *, integer *), 
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

#line 194 "cpbstf.f"
    /* Parameter adjustments */
#line 194 "cpbstf.f"
    ab_dim1 = *ldab;
#line 194 "cpbstf.f"
    ab_offset = 1 + ab_dim1;
#line 194 "cpbstf.f"
    ab -= ab_offset;
#line 194 "cpbstf.f"

#line 194 "cpbstf.f"
    /* Function Body */
#line 194 "cpbstf.f"
    *info = 0;
#line 195 "cpbstf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 196 "cpbstf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 197 "cpbstf.f"
	*info = -1;
#line 198 "cpbstf.f"
    } else if (*n < 0) {
#line 199 "cpbstf.f"
	*info = -2;
#line 200 "cpbstf.f"
    } else if (*kd < 0) {
#line 201 "cpbstf.f"
	*info = -3;
#line 202 "cpbstf.f"
    } else if (*ldab < *kd + 1) {
#line 203 "cpbstf.f"
	*info = -5;
#line 204 "cpbstf.f"
    }
#line 205 "cpbstf.f"
    if (*info != 0) {
#line 206 "cpbstf.f"
	i__1 = -(*info);
#line 206 "cpbstf.f"
	xerbla_("CPBSTF", &i__1, (ftnlen)6);
#line 207 "cpbstf.f"
	return 0;
#line 208 "cpbstf.f"
    }

/*     Quick return if possible */

#line 212 "cpbstf.f"
    if (*n == 0) {
#line 212 "cpbstf.f"
	return 0;
#line 212 "cpbstf.f"
    }

/* Computing MAX */
#line 215 "cpbstf.f"
    i__1 = 1, i__2 = *ldab - 1;
#line 215 "cpbstf.f"
    kld = max(i__1,i__2);

/*     Set the splitting point m. */

#line 219 "cpbstf.f"
    m = (*n + *kd) / 2;

#line 221 "cpbstf.f"
    if (upper) {

/*        Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m). */

#line 225 "cpbstf.f"
	i__1 = m + 1;
#line 225 "cpbstf.f"
	for (j = *n; j >= i__1; --j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 229 "cpbstf.f"
	    i__2 = *kd + 1 + j * ab_dim1;
#line 229 "cpbstf.f"
	    ajj = ab[i__2].r;
#line 230 "cpbstf.f"
	    if (ajj <= 0.) {
#line 231 "cpbstf.f"
		i__2 = *kd + 1 + j * ab_dim1;
#line 231 "cpbstf.f"
		ab[i__2].r = ajj, ab[i__2].i = 0.;
#line 232 "cpbstf.f"
		goto L50;
#line 233 "cpbstf.f"
	    }
#line 234 "cpbstf.f"
	    ajj = sqrt(ajj);
#line 235 "cpbstf.f"
	    i__2 = *kd + 1 + j * ab_dim1;
#line 235 "cpbstf.f"
	    ab[i__2].r = ajj, ab[i__2].i = 0.;
/* Computing MIN */
#line 236 "cpbstf.f"
	    i__2 = j - 1;
#line 236 "cpbstf.f"
	    km = min(i__2,*kd);

/*           Compute elements j-km:j-1 of the j-th column and update the */
/*           the leading submatrix within the band. */

#line 241 "cpbstf.f"
	    d__1 = 1. / ajj;
#line 241 "cpbstf.f"
	    csscal_(&km, &d__1, &ab[*kd + 1 - km + j * ab_dim1], &c__1);
#line 242 "cpbstf.f"
	    cher_("Upper", &km, &c_b9, &ab[*kd + 1 - km + j * ab_dim1], &c__1,
		     &ab[*kd + 1 + (j - km) * ab_dim1], &kld, (ftnlen)5);
#line 244 "cpbstf.f"
/* L10: */
#line 244 "cpbstf.f"
	}

/*        Factorize the updated submatrix A(1:m,1:m) as U**H*U. */

#line 248 "cpbstf.f"
	i__1 = m;
#line 248 "cpbstf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 252 "cpbstf.f"
	    i__2 = *kd + 1 + j * ab_dim1;
#line 252 "cpbstf.f"
	    ajj = ab[i__2].r;
#line 253 "cpbstf.f"
	    if (ajj <= 0.) {
#line 254 "cpbstf.f"
		i__2 = *kd + 1 + j * ab_dim1;
#line 254 "cpbstf.f"
		ab[i__2].r = ajj, ab[i__2].i = 0.;
#line 255 "cpbstf.f"
		goto L50;
#line 256 "cpbstf.f"
	    }
#line 257 "cpbstf.f"
	    ajj = sqrt(ajj);
#line 258 "cpbstf.f"
	    i__2 = *kd + 1 + j * ab_dim1;
#line 258 "cpbstf.f"
	    ab[i__2].r = ajj, ab[i__2].i = 0.;
/* Computing MIN */
#line 259 "cpbstf.f"
	    i__2 = *kd, i__3 = m - j;
#line 259 "cpbstf.f"
	    km = min(i__2,i__3);

/*           Compute elements j+1:j+km of the j-th row and update the */
/*           trailing submatrix within the band. */

#line 264 "cpbstf.f"
	    if (km > 0) {
#line 265 "cpbstf.f"
		d__1 = 1. / ajj;
#line 265 "cpbstf.f"
		csscal_(&km, &d__1, &ab[*kd + (j + 1) * ab_dim1], &kld);
#line 266 "cpbstf.f"
		clacgv_(&km, &ab[*kd + (j + 1) * ab_dim1], &kld);
#line 267 "cpbstf.f"
		cher_("Upper", &km, &c_b9, &ab[*kd + (j + 1) * ab_dim1], &kld,
			 &ab[*kd + 1 + (j + 1) * ab_dim1], &kld, (ftnlen)5);
#line 269 "cpbstf.f"
		clacgv_(&km, &ab[*kd + (j + 1) * ab_dim1], &kld);
#line 270 "cpbstf.f"
	    }
#line 271 "cpbstf.f"
/* L20: */
#line 271 "cpbstf.f"
	}
#line 272 "cpbstf.f"
    } else {

/*        Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m). */

#line 276 "cpbstf.f"
	i__1 = m + 1;
#line 276 "cpbstf.f"
	for (j = *n; j >= i__1; --j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 280 "cpbstf.f"
	    i__2 = j * ab_dim1 + 1;
#line 280 "cpbstf.f"
	    ajj = ab[i__2].r;
#line 281 "cpbstf.f"
	    if (ajj <= 0.) {
#line 282 "cpbstf.f"
		i__2 = j * ab_dim1 + 1;
#line 282 "cpbstf.f"
		ab[i__2].r = ajj, ab[i__2].i = 0.;
#line 283 "cpbstf.f"
		goto L50;
#line 284 "cpbstf.f"
	    }
#line 285 "cpbstf.f"
	    ajj = sqrt(ajj);
#line 286 "cpbstf.f"
	    i__2 = j * ab_dim1 + 1;
#line 286 "cpbstf.f"
	    ab[i__2].r = ajj, ab[i__2].i = 0.;
/* Computing MIN */
#line 287 "cpbstf.f"
	    i__2 = j - 1;
#line 287 "cpbstf.f"
	    km = min(i__2,*kd);

/*           Compute elements j-km:j-1 of the j-th row and update the */
/*           trailing submatrix within the band. */

#line 292 "cpbstf.f"
	    d__1 = 1. / ajj;
#line 292 "cpbstf.f"
	    csscal_(&km, &d__1, &ab[km + 1 + (j - km) * ab_dim1], &kld);
#line 293 "cpbstf.f"
	    clacgv_(&km, &ab[km + 1 + (j - km) * ab_dim1], &kld);
#line 294 "cpbstf.f"
	    cher_("Lower", &km, &c_b9, &ab[km + 1 + (j - km) * ab_dim1], &kld,
		     &ab[(j - km) * ab_dim1 + 1], &kld, (ftnlen)5);
#line 296 "cpbstf.f"
	    clacgv_(&km, &ab[km + 1 + (j - km) * ab_dim1], &kld);
#line 297 "cpbstf.f"
/* L30: */
#line 297 "cpbstf.f"
	}

/*        Factorize the updated submatrix A(1:m,1:m) as U**H*U. */

#line 301 "cpbstf.f"
	i__1 = m;
#line 301 "cpbstf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 305 "cpbstf.f"
	    i__2 = j * ab_dim1 + 1;
#line 305 "cpbstf.f"
	    ajj = ab[i__2].r;
#line 306 "cpbstf.f"
	    if (ajj <= 0.) {
#line 307 "cpbstf.f"
		i__2 = j * ab_dim1 + 1;
#line 307 "cpbstf.f"
		ab[i__2].r = ajj, ab[i__2].i = 0.;
#line 308 "cpbstf.f"
		goto L50;
#line 309 "cpbstf.f"
	    }
#line 310 "cpbstf.f"
	    ajj = sqrt(ajj);
#line 311 "cpbstf.f"
	    i__2 = j * ab_dim1 + 1;
#line 311 "cpbstf.f"
	    ab[i__2].r = ajj, ab[i__2].i = 0.;
/* Computing MIN */
#line 312 "cpbstf.f"
	    i__2 = *kd, i__3 = m - j;
#line 312 "cpbstf.f"
	    km = min(i__2,i__3);

/*           Compute elements j+1:j+km of the j-th column and update the */
/*           trailing submatrix within the band. */

#line 317 "cpbstf.f"
	    if (km > 0) {
#line 318 "cpbstf.f"
		d__1 = 1. / ajj;
#line 318 "cpbstf.f"
		csscal_(&km, &d__1, &ab[j * ab_dim1 + 2], &c__1);
#line 319 "cpbstf.f"
		cher_("Lower", &km, &c_b9, &ab[j * ab_dim1 + 2], &c__1, &ab[(
			j + 1) * ab_dim1 + 1], &kld, (ftnlen)5);
#line 321 "cpbstf.f"
	    }
#line 322 "cpbstf.f"
/* L40: */
#line 322 "cpbstf.f"
	}
#line 323 "cpbstf.f"
    }
#line 324 "cpbstf.f"
    return 0;

#line 326 "cpbstf.f"
L50:
#line 327 "cpbstf.f"
    *info = j;
#line 328 "cpbstf.f"
    return 0;

/*     End of CPBSTF */

} /* cpbstf_ */

