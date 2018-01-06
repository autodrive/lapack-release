#line 1 "spbstf.f"
/* spbstf.f -- translated by f2c (version 20100827).
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

#line 1 "spbstf.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = -1.;

/* > \brief \b SPBSTF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPBSTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spbstf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spbstf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spbstf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPBSTF( UPLO, N, KD, AB, LDAB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPBSTF computes a split Cholesky factorization of a real */
/* > symmetric positive definite band matrix A. */
/* > */
/* > This routine is designed to be used in conjunction with SSBGST. */
/* > */
/* > The factorization has the form  A = S**T*S  where S is a band matrix */
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
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first kd+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, if INFO = 0, the factor S from the split Cholesky */
/* >          factorization A = S**T*S. See Further Details. */
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

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

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
/* >   *    *   a13  a24  a35  a46  a57   *    *   s13  s24  s53  s64  s75 */
/* >   *   a12  a23  a34  a45  a56  a67   *   s12  s23  s34  s54  s65  s76 */
/* >  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55  s66  s77 */
/* > */
/* >  If UPLO = 'L', the array AB holds: */
/* > */
/* >  on entry:                          on exit: */
/* > */
/* >  a11  a22  a33  a44  a55  a66  a77  s11  s22  s33  s44  s55  s66  s77 */
/* >  a21  a32  a43  a54  a65  a76   *   s12  s23  s34  s54  s65  s76   * */
/* >  a31  a42  a53  a64  a64   *    *   s13  s24  s53  s64  s75   *    * */
/* > */
/* >  Array elements marked * are not used by the routine. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int spbstf_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info, ftnlen uplo_len)
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
    extern /* Subroutine */ int ssyr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

#line 193 "spbstf.f"
    /* Parameter adjustments */
#line 193 "spbstf.f"
    ab_dim1 = *ldab;
#line 193 "spbstf.f"
    ab_offset = 1 + ab_dim1;
#line 193 "spbstf.f"
    ab -= ab_offset;
#line 193 "spbstf.f"

#line 193 "spbstf.f"
    /* Function Body */
#line 193 "spbstf.f"
    *info = 0;
#line 194 "spbstf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 195 "spbstf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 196 "spbstf.f"
	*info = -1;
#line 197 "spbstf.f"
    } else if (*n < 0) {
#line 198 "spbstf.f"
	*info = -2;
#line 199 "spbstf.f"
    } else if (*kd < 0) {
#line 200 "spbstf.f"
	*info = -3;
#line 201 "spbstf.f"
    } else if (*ldab < *kd + 1) {
#line 202 "spbstf.f"
	*info = -5;
#line 203 "spbstf.f"
    }
#line 204 "spbstf.f"
    if (*info != 0) {
#line 205 "spbstf.f"
	i__1 = -(*info);
#line 205 "spbstf.f"
	xerbla_("SPBSTF", &i__1, (ftnlen)6);
#line 206 "spbstf.f"
	return 0;
#line 207 "spbstf.f"
    }

/*     Quick return if possible */

#line 211 "spbstf.f"
    if (*n == 0) {
#line 211 "spbstf.f"
	return 0;
#line 211 "spbstf.f"
    }

/* Computing MAX */
#line 214 "spbstf.f"
    i__1 = 1, i__2 = *ldab - 1;
#line 214 "spbstf.f"
    kld = max(i__1,i__2);

/*     Set the splitting point m. */

#line 218 "spbstf.f"
    m = (*n + *kd) / 2;

#line 220 "spbstf.f"
    if (upper) {

/*        Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m). */

#line 224 "spbstf.f"
	i__1 = m + 1;
#line 224 "spbstf.f"
	for (j = *n; j >= i__1; --j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 228 "spbstf.f"
	    ajj = ab[*kd + 1 + j * ab_dim1];
#line 229 "spbstf.f"
	    if (ajj <= 0.) {
#line 229 "spbstf.f"
		goto L50;
#line 229 "spbstf.f"
	    }
#line 231 "spbstf.f"
	    ajj = sqrt(ajj);
#line 232 "spbstf.f"
	    ab[*kd + 1 + j * ab_dim1] = ajj;
/* Computing MIN */
#line 233 "spbstf.f"
	    i__2 = j - 1;
#line 233 "spbstf.f"
	    km = min(i__2,*kd);

/*           Compute elements j-km:j-1 of the j-th column and update the */
/*           the leading submatrix within the band. */

#line 238 "spbstf.f"
	    d__1 = 1. / ajj;
#line 238 "spbstf.f"
	    sscal_(&km, &d__1, &ab[*kd + 1 - km + j * ab_dim1], &c__1);
#line 239 "spbstf.f"
	    ssyr_("Upper", &km, &c_b9, &ab[*kd + 1 - km + j * ab_dim1], &c__1,
		     &ab[*kd + 1 + (j - km) * ab_dim1], &kld, (ftnlen)5);
#line 241 "spbstf.f"
/* L10: */
#line 241 "spbstf.f"
	}

/*        Factorize the updated submatrix A(1:m,1:m) as U**T*U. */

#line 245 "spbstf.f"
	i__1 = m;
#line 245 "spbstf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 249 "spbstf.f"
	    ajj = ab[*kd + 1 + j * ab_dim1];
#line 250 "spbstf.f"
	    if (ajj <= 0.) {
#line 250 "spbstf.f"
		goto L50;
#line 250 "spbstf.f"
	    }
#line 252 "spbstf.f"
	    ajj = sqrt(ajj);
#line 253 "spbstf.f"
	    ab[*kd + 1 + j * ab_dim1] = ajj;
/* Computing MIN */
#line 254 "spbstf.f"
	    i__2 = *kd, i__3 = m - j;
#line 254 "spbstf.f"
	    km = min(i__2,i__3);

/*           Compute elements j+1:j+km of the j-th row and update the */
/*           trailing submatrix within the band. */

#line 259 "spbstf.f"
	    if (km > 0) {
#line 260 "spbstf.f"
		d__1 = 1. / ajj;
#line 260 "spbstf.f"
		sscal_(&km, &d__1, &ab[*kd + (j + 1) * ab_dim1], &kld);
#line 261 "spbstf.f"
		ssyr_("Upper", &km, &c_b9, &ab[*kd + (j + 1) * ab_dim1], &kld,
			 &ab[*kd + 1 + (j + 1) * ab_dim1], &kld, (ftnlen)5);
#line 263 "spbstf.f"
	    }
#line 264 "spbstf.f"
/* L20: */
#line 264 "spbstf.f"
	}
#line 265 "spbstf.f"
    } else {

/*        Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m). */

#line 269 "spbstf.f"
	i__1 = m + 1;
#line 269 "spbstf.f"
	for (j = *n; j >= i__1; --j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 273 "spbstf.f"
	    ajj = ab[j * ab_dim1 + 1];
#line 274 "spbstf.f"
	    if (ajj <= 0.) {
#line 274 "spbstf.f"
		goto L50;
#line 274 "spbstf.f"
	    }
#line 276 "spbstf.f"
	    ajj = sqrt(ajj);
#line 277 "spbstf.f"
	    ab[j * ab_dim1 + 1] = ajj;
/* Computing MIN */
#line 278 "spbstf.f"
	    i__2 = j - 1;
#line 278 "spbstf.f"
	    km = min(i__2,*kd);

/*           Compute elements j-km:j-1 of the j-th row and update the */
/*           trailing submatrix within the band. */

#line 283 "spbstf.f"
	    d__1 = 1. / ajj;
#line 283 "spbstf.f"
	    sscal_(&km, &d__1, &ab[km + 1 + (j - km) * ab_dim1], &kld);
#line 284 "spbstf.f"
	    ssyr_("Lower", &km, &c_b9, &ab[km + 1 + (j - km) * ab_dim1], &kld,
		     &ab[(j - km) * ab_dim1 + 1], &kld, (ftnlen)5);
#line 286 "spbstf.f"
/* L30: */
#line 286 "spbstf.f"
	}

/*        Factorize the updated submatrix A(1:m,1:m) as U**T*U. */

#line 290 "spbstf.f"
	i__1 = m;
#line 290 "spbstf.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute s(j,j) and test for non-positive-definiteness. */

#line 294 "spbstf.f"
	    ajj = ab[j * ab_dim1 + 1];
#line 295 "spbstf.f"
	    if (ajj <= 0.) {
#line 295 "spbstf.f"
		goto L50;
#line 295 "spbstf.f"
	    }
#line 297 "spbstf.f"
	    ajj = sqrt(ajj);
#line 298 "spbstf.f"
	    ab[j * ab_dim1 + 1] = ajj;
/* Computing MIN */
#line 299 "spbstf.f"
	    i__2 = *kd, i__3 = m - j;
#line 299 "spbstf.f"
	    km = min(i__2,i__3);

/*           Compute elements j+1:j+km of the j-th column and update the */
/*           trailing submatrix within the band. */

#line 304 "spbstf.f"
	    if (km > 0) {
#line 305 "spbstf.f"
		d__1 = 1. / ajj;
#line 305 "spbstf.f"
		sscal_(&km, &d__1, &ab[j * ab_dim1 + 2], &c__1);
#line 306 "spbstf.f"
		ssyr_("Lower", &km, &c_b9, &ab[j * ab_dim1 + 2], &c__1, &ab[(
			j + 1) * ab_dim1 + 1], &kld, (ftnlen)5);
#line 308 "spbstf.f"
	    }
#line 309 "spbstf.f"
/* L40: */
#line 309 "spbstf.f"
	}
#line 310 "spbstf.f"
    }
#line 311 "spbstf.f"
    return 0;

#line 313 "spbstf.f"
L50:
#line 314 "spbstf.f"
    *info = j;
#line 315 "spbstf.f"
    return 0;

/*     End of SPBSTF */

} /* spbstf_ */

