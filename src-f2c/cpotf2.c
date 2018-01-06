#line 1 "cpotf2.f"
/* cpotf2.f -- translated by f2c (version 20100827).
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

#line 1 "cpotf2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CPOTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite matrix (u
nblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPOTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpotf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpotf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpotf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPOTF2( UPLO, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPOTF2 computes the Cholesky factorization of a complex Hermitian */
/* > positive definite matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**H * U ,  if UPLO = 'U', or */
/* >    A = L  * L**H,  if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
/* > */
/* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          Hermitian matrix A is stored. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
/* >          n by n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n by n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization A = U**H *U  or A = L*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
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

/* > \ingroup complexPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpotf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j;
    static doublereal ajj;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , csscal_(integer *, doublereal *, doublecomplex *, integer *), 
	    xerbla_(char *, integer *, ftnlen);
    extern logical sisnan_(doublereal *);


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

#line 153 "cpotf2.f"
    /* Parameter adjustments */
#line 153 "cpotf2.f"
    a_dim1 = *lda;
#line 153 "cpotf2.f"
    a_offset = 1 + a_dim1;
#line 153 "cpotf2.f"
    a -= a_offset;
#line 153 "cpotf2.f"

#line 153 "cpotf2.f"
    /* Function Body */
#line 153 "cpotf2.f"
    *info = 0;
#line 154 "cpotf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 155 "cpotf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 156 "cpotf2.f"
	*info = -1;
#line 157 "cpotf2.f"
    } else if (*n < 0) {
#line 158 "cpotf2.f"
	*info = -2;
#line 159 "cpotf2.f"
    } else if (*lda < max(1,*n)) {
#line 160 "cpotf2.f"
	*info = -4;
#line 161 "cpotf2.f"
    }
#line 162 "cpotf2.f"
    if (*info != 0) {
#line 163 "cpotf2.f"
	i__1 = -(*info);
#line 163 "cpotf2.f"
	xerbla_("CPOTF2", &i__1, (ftnlen)6);
#line 164 "cpotf2.f"
	return 0;
#line 165 "cpotf2.f"
    }

/*     Quick return if possible */

#line 169 "cpotf2.f"
    if (*n == 0) {
#line 169 "cpotf2.f"
	return 0;
#line 169 "cpotf2.f"
    }

#line 172 "cpotf2.f"
    if (upper) {

/*        Compute the Cholesky factorization A = U**H *U. */

#line 176 "cpotf2.f"
	i__1 = *n;
#line 176 "cpotf2.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness. */

#line 180 "cpotf2.f"
	    i__2 = j + j * a_dim1;
#line 180 "cpotf2.f"
	    d__1 = a[i__2].r;
#line 180 "cpotf2.f"
	    i__3 = j - 1;
#line 180 "cpotf2.f"
	    cdotc_(&z__2, &i__3, &a[j * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1]
		    , &c__1);
#line 180 "cpotf2.f"
	    z__1.r = d__1 - z__2.r, z__1.i = -z__2.i;
#line 180 "cpotf2.f"
	    ajj = z__1.r;
#line 182 "cpotf2.f"
	    if (ajj <= 0. || sisnan_(&ajj)) {
#line 183 "cpotf2.f"
		i__2 = j + j * a_dim1;
#line 183 "cpotf2.f"
		a[i__2].r = ajj, a[i__2].i = 0.;
#line 184 "cpotf2.f"
		goto L30;
#line 185 "cpotf2.f"
	    }
#line 186 "cpotf2.f"
	    ajj = sqrt(ajj);
#line 187 "cpotf2.f"
	    i__2 = j + j * a_dim1;
#line 187 "cpotf2.f"
	    a[i__2].r = ajj, a[i__2].i = 0.;

/*           Compute elements J+1:N of row J. */

#line 191 "cpotf2.f"
	    if (j < *n) {
#line 192 "cpotf2.f"
		i__2 = j - 1;
#line 192 "cpotf2.f"
		clacgv_(&i__2, &a[j * a_dim1 + 1], &c__1);
#line 193 "cpotf2.f"
		i__2 = j - 1;
#line 193 "cpotf2.f"
		i__3 = *n - j;
#line 193 "cpotf2.f"
		z__1.r = -1., z__1.i = -0.;
#line 193 "cpotf2.f"
		cgemv_("Transpose", &i__2, &i__3, &z__1, &a[(j + 1) * a_dim1 
			+ 1], lda, &a[j * a_dim1 + 1], &c__1, &c_b1, &a[j + (
			j + 1) * a_dim1], lda, (ftnlen)9);
#line 195 "cpotf2.f"
		i__2 = j - 1;
#line 195 "cpotf2.f"
		clacgv_(&i__2, &a[j * a_dim1 + 1], &c__1);
#line 196 "cpotf2.f"
		i__2 = *n - j;
#line 196 "cpotf2.f"
		d__1 = 1. / ajj;
#line 196 "cpotf2.f"
		csscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
#line 197 "cpotf2.f"
	    }
#line 198 "cpotf2.f"
/* L10: */
#line 198 "cpotf2.f"
	}
#line 199 "cpotf2.f"
    } else {

/*        Compute the Cholesky factorization A = L*L**H. */

#line 203 "cpotf2.f"
	i__1 = *n;
#line 203 "cpotf2.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

#line 207 "cpotf2.f"
	    i__2 = j + j * a_dim1;
#line 207 "cpotf2.f"
	    d__1 = a[i__2].r;
#line 207 "cpotf2.f"
	    i__3 = j - 1;
#line 207 "cpotf2.f"
	    cdotc_(&z__2, &i__3, &a[j + a_dim1], lda, &a[j + a_dim1], lda);
#line 207 "cpotf2.f"
	    z__1.r = d__1 - z__2.r, z__1.i = -z__2.i;
#line 207 "cpotf2.f"
	    ajj = z__1.r;
#line 209 "cpotf2.f"
	    if (ajj <= 0. || sisnan_(&ajj)) {
#line 210 "cpotf2.f"
		i__2 = j + j * a_dim1;
#line 210 "cpotf2.f"
		a[i__2].r = ajj, a[i__2].i = 0.;
#line 211 "cpotf2.f"
		goto L30;
#line 212 "cpotf2.f"
	    }
#line 213 "cpotf2.f"
	    ajj = sqrt(ajj);
#line 214 "cpotf2.f"
	    i__2 = j + j * a_dim1;
#line 214 "cpotf2.f"
	    a[i__2].r = ajj, a[i__2].i = 0.;

/*           Compute elements J+1:N of column J. */

#line 218 "cpotf2.f"
	    if (j < *n) {
#line 219 "cpotf2.f"
		i__2 = j - 1;
#line 219 "cpotf2.f"
		clacgv_(&i__2, &a[j + a_dim1], lda);
#line 220 "cpotf2.f"
		i__2 = *n - j;
#line 220 "cpotf2.f"
		i__3 = j - 1;
#line 220 "cpotf2.f"
		z__1.r = -1., z__1.i = -0.;
#line 220 "cpotf2.f"
		cgemv_("No transpose", &i__2, &i__3, &z__1, &a[j + 1 + a_dim1]
			, lda, &a[j + a_dim1], lda, &c_b1, &a[j + 1 + j * 
			a_dim1], &c__1, (ftnlen)12);
#line 222 "cpotf2.f"
		i__2 = j - 1;
#line 222 "cpotf2.f"
		clacgv_(&i__2, &a[j + a_dim1], lda);
#line 223 "cpotf2.f"
		i__2 = *n - j;
#line 223 "cpotf2.f"
		d__1 = 1. / ajj;
#line 223 "cpotf2.f"
		csscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
#line 224 "cpotf2.f"
	    }
#line 225 "cpotf2.f"
/* L20: */
#line 225 "cpotf2.f"
	}
#line 226 "cpotf2.f"
    }
#line 227 "cpotf2.f"
    goto L40;

#line 229 "cpotf2.f"
L30:
#line 230 "cpotf2.f"
    *info = j;

#line 232 "cpotf2.f"
L40:
#line 233 "cpotf2.f"
    return 0;

/*     End of CPOTF2 */

} /* cpotf2_ */

