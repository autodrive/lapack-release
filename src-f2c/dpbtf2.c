#line 1 "dpbtf2.f"
/* dpbtf2.f -- translated by f2c (version 20100827).
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

#line 1 "dpbtf2.f"
/* Table of constant values */

static doublereal c_b8 = -1.;
static integer c__1 = 1;

/* > \brief \b DPBTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite band matr
ix (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPBTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbtf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbtf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbtf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPBTF2 computes the Cholesky factorization of a real symmetric */
/* > positive definite band matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**T * U ,  if UPLO = 'U', or */
/* >    A = L  * L**T,  if UPLO = 'L', */
/* > where U is an upper triangular matrix, U**T is the transpose of U, and */
/* > L is lower triangular. */
/* > */
/* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored: */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
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
/* >          The number of super-diagonals of the matrix A if UPLO = 'U', */
/* >          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, if INFO = 0, the triangular factor U or L from the */
/* >          Cholesky factorization A = U**T*U or A = L*L**T of the band */
/* >          matrix A, in the same storage format as A. */
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
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* >          > 0: if INFO = k, the leading minor of order k is not */
/* >               positive definite, and the factorization could not be */
/* >               completed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The band storage scheme is illustrated by the following example, when */
/* >  N = 6, KD = 2, and UPLO = 'U': */
/* > */
/* >  On entry:                       On exit: */
/* > */
/* >      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46 */
/* >      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56 */
/* >     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66 */
/* > */
/* >  Similarly, if UPLO = 'L' the format of A is as follows: */
/* > */
/* >  On entry:                       On exit: */
/* > */
/* >     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66 */
/* >     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   * */
/* >     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    * */
/* > */
/* >  Array elements marked * are not used by the routine. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dpbtf2_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, kn;
    static doublereal ajj;
    static integer kld;
    extern /* Subroutine */ int dsyr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 183 "dpbtf2.f"
    /* Parameter adjustments */
#line 183 "dpbtf2.f"
    ab_dim1 = *ldab;
#line 183 "dpbtf2.f"
    ab_offset = 1 + ab_dim1;
#line 183 "dpbtf2.f"
    ab -= ab_offset;
#line 183 "dpbtf2.f"

#line 183 "dpbtf2.f"
    /* Function Body */
#line 183 "dpbtf2.f"
    *info = 0;
#line 184 "dpbtf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 185 "dpbtf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 186 "dpbtf2.f"
	*info = -1;
#line 187 "dpbtf2.f"
    } else if (*n < 0) {
#line 188 "dpbtf2.f"
	*info = -2;
#line 189 "dpbtf2.f"
    } else if (*kd < 0) {
#line 190 "dpbtf2.f"
	*info = -3;
#line 191 "dpbtf2.f"
    } else if (*ldab < *kd + 1) {
#line 192 "dpbtf2.f"
	*info = -5;
#line 193 "dpbtf2.f"
    }
#line 194 "dpbtf2.f"
    if (*info != 0) {
#line 195 "dpbtf2.f"
	i__1 = -(*info);
#line 195 "dpbtf2.f"
	xerbla_("DPBTF2", &i__1, (ftnlen)6);
#line 196 "dpbtf2.f"
	return 0;
#line 197 "dpbtf2.f"
    }

/*     Quick return if possible */

#line 201 "dpbtf2.f"
    if (*n == 0) {
#line 201 "dpbtf2.f"
	return 0;
#line 201 "dpbtf2.f"
    }

/* Computing MAX */
#line 204 "dpbtf2.f"
    i__1 = 1, i__2 = *ldab - 1;
#line 204 "dpbtf2.f"
    kld = max(i__1,i__2);

#line 206 "dpbtf2.f"
    if (upper) {

/*        Compute the Cholesky factorization A = U**T*U. */

#line 210 "dpbtf2.f"
	i__1 = *n;
#line 210 "dpbtf2.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness. */

#line 214 "dpbtf2.f"
	    ajj = ab[*kd + 1 + j * ab_dim1];
#line 215 "dpbtf2.f"
	    if (ajj <= 0.) {
#line 215 "dpbtf2.f"
		goto L30;
#line 215 "dpbtf2.f"
	    }
#line 217 "dpbtf2.f"
	    ajj = sqrt(ajj);
#line 218 "dpbtf2.f"
	    ab[*kd + 1 + j * ab_dim1] = ajj;

/*           Compute elements J+1:J+KN of row J and update the */
/*           trailing submatrix within the band. */

/* Computing MIN */
#line 223 "dpbtf2.f"
	    i__2 = *kd, i__3 = *n - j;
#line 223 "dpbtf2.f"
	    kn = min(i__2,i__3);
#line 224 "dpbtf2.f"
	    if (kn > 0) {
#line 225 "dpbtf2.f"
		d__1 = 1. / ajj;
#line 225 "dpbtf2.f"
		dscal_(&kn, &d__1, &ab[*kd + (j + 1) * ab_dim1], &kld);
#line 226 "dpbtf2.f"
		dsyr_("Upper", &kn, &c_b8, &ab[*kd + (j + 1) * ab_dim1], &kld,
			 &ab[*kd + 1 + (j + 1) * ab_dim1], &kld, (ftnlen)5);
#line 228 "dpbtf2.f"
	    }
#line 229 "dpbtf2.f"
/* L10: */
#line 229 "dpbtf2.f"
	}
#line 230 "dpbtf2.f"
    } else {

/*        Compute the Cholesky factorization A = L*L**T. */

#line 234 "dpbtf2.f"
	i__1 = *n;
#line 234 "dpbtf2.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

#line 238 "dpbtf2.f"
	    ajj = ab[j * ab_dim1 + 1];
#line 239 "dpbtf2.f"
	    if (ajj <= 0.) {
#line 239 "dpbtf2.f"
		goto L30;
#line 239 "dpbtf2.f"
	    }
#line 241 "dpbtf2.f"
	    ajj = sqrt(ajj);
#line 242 "dpbtf2.f"
	    ab[j * ab_dim1 + 1] = ajj;

/*           Compute elements J+1:J+KN of column J and update the */
/*           trailing submatrix within the band. */

/* Computing MIN */
#line 247 "dpbtf2.f"
	    i__2 = *kd, i__3 = *n - j;
#line 247 "dpbtf2.f"
	    kn = min(i__2,i__3);
#line 248 "dpbtf2.f"
	    if (kn > 0) {
#line 249 "dpbtf2.f"
		d__1 = 1. / ajj;
#line 249 "dpbtf2.f"
		dscal_(&kn, &d__1, &ab[j * ab_dim1 + 2], &c__1);
#line 250 "dpbtf2.f"
		dsyr_("Lower", &kn, &c_b8, &ab[j * ab_dim1 + 2], &c__1, &ab[(
			j + 1) * ab_dim1 + 1], &kld, (ftnlen)5);
#line 252 "dpbtf2.f"
	    }
#line 253 "dpbtf2.f"
/* L20: */
#line 253 "dpbtf2.f"
	}
#line 254 "dpbtf2.f"
    }
#line 255 "dpbtf2.f"
    return 0;

#line 257 "dpbtf2.f"
L30:
#line 258 "dpbtf2.f"
    *info = j;
#line 259 "dpbtf2.f"
    return 0;

/*     End of DPBTF2 */

} /* dpbtf2_ */

