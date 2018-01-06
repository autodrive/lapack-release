#line 1 "spotf2.f"
/* spotf2.f -- translated by f2c (version 20100827).
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

#line 1 "spotf2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b10 = -1.;
static doublereal c_b12 = 1.;

/* > \brief \b SPOTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite matrix (u
nblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPOTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spotf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spotf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spotf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPOTF2( UPLO, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPOTF2 computes the Cholesky factorization of a real symmetric */
/* > positive definite matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**T * U ,  if UPLO = 'U', or */
/* >    A = L  * L**T,  if UPLO = 'L', */
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
/* >          symmetric matrix A is stored. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          n by n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n by n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization A = U**T *U  or A = L*L**T. */
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

/* > \date December 2016 */

/* > \ingroup realPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int spotf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j;
    static doublereal ajj;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern logical sisnan_(doublereal *);


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

#line 151 "spotf2.f"
    /* Parameter adjustments */
#line 151 "spotf2.f"
    a_dim1 = *lda;
#line 151 "spotf2.f"
    a_offset = 1 + a_dim1;
#line 151 "spotf2.f"
    a -= a_offset;
#line 151 "spotf2.f"

#line 151 "spotf2.f"
    /* Function Body */
#line 151 "spotf2.f"
    *info = 0;
#line 152 "spotf2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 153 "spotf2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 154 "spotf2.f"
	*info = -1;
#line 155 "spotf2.f"
    } else if (*n < 0) {
#line 156 "spotf2.f"
	*info = -2;
#line 157 "spotf2.f"
    } else if (*lda < max(1,*n)) {
#line 158 "spotf2.f"
	*info = -4;
#line 159 "spotf2.f"
    }
#line 160 "spotf2.f"
    if (*info != 0) {
#line 161 "spotf2.f"
	i__1 = -(*info);
#line 161 "spotf2.f"
	xerbla_("SPOTF2", &i__1, (ftnlen)6);
#line 162 "spotf2.f"
	return 0;
#line 163 "spotf2.f"
    }

/*     Quick return if possible */

#line 167 "spotf2.f"
    if (*n == 0) {
#line 167 "spotf2.f"
	return 0;
#line 167 "spotf2.f"
    }

#line 170 "spotf2.f"
    if (upper) {

/*        Compute the Cholesky factorization A = U**T *U. */

#line 174 "spotf2.f"
	i__1 = *n;
#line 174 "spotf2.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness. */

#line 178 "spotf2.f"
	    i__2 = j - 1;
#line 178 "spotf2.f"
	    ajj = a[j + j * a_dim1] - sdot_(&i__2, &a[j * a_dim1 + 1], &c__1, 
		    &a[j * a_dim1 + 1], &c__1);
#line 179 "spotf2.f"
	    if (ajj <= 0. || sisnan_(&ajj)) {
#line 180 "spotf2.f"
		a[j + j * a_dim1] = ajj;
#line 181 "spotf2.f"
		goto L30;
#line 182 "spotf2.f"
	    }
#line 183 "spotf2.f"
	    ajj = sqrt(ajj);
#line 184 "spotf2.f"
	    a[j + j * a_dim1] = ajj;

/*           Compute elements J+1:N of row J. */

#line 188 "spotf2.f"
	    if (j < *n) {
#line 189 "spotf2.f"
		i__2 = j - 1;
#line 189 "spotf2.f"
		i__3 = *n - j;
#line 189 "spotf2.f"
		sgemv_("Transpose", &i__2, &i__3, &c_b10, &a[(j + 1) * a_dim1 
			+ 1], lda, &a[j * a_dim1 + 1], &c__1, &c_b12, &a[j + (
			j + 1) * a_dim1], lda, (ftnlen)9);
#line 191 "spotf2.f"
		i__2 = *n - j;
#line 191 "spotf2.f"
		d__1 = 1. / ajj;
#line 191 "spotf2.f"
		sscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
#line 192 "spotf2.f"
	    }
#line 193 "spotf2.f"
/* L10: */
#line 193 "spotf2.f"
	}
#line 194 "spotf2.f"
    } else {

/*        Compute the Cholesky factorization A = L*L**T. */

#line 198 "spotf2.f"
	i__1 = *n;
#line 198 "spotf2.f"
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

#line 202 "spotf2.f"
	    i__2 = j - 1;
#line 202 "spotf2.f"
	    ajj = a[j + j * a_dim1] - sdot_(&i__2, &a[j + a_dim1], lda, &a[j 
		    + a_dim1], lda);
#line 204 "spotf2.f"
	    if (ajj <= 0. || sisnan_(&ajj)) {
#line 205 "spotf2.f"
		a[j + j * a_dim1] = ajj;
#line 206 "spotf2.f"
		goto L30;
#line 207 "spotf2.f"
	    }
#line 208 "spotf2.f"
	    ajj = sqrt(ajj);
#line 209 "spotf2.f"
	    a[j + j * a_dim1] = ajj;

/*           Compute elements J+1:N of column J. */

#line 213 "spotf2.f"
	    if (j < *n) {
#line 214 "spotf2.f"
		i__2 = *n - j;
#line 214 "spotf2.f"
		i__3 = j - 1;
#line 214 "spotf2.f"
		sgemv_("No transpose", &i__2, &i__3, &c_b10, &a[j + 1 + 
			a_dim1], lda, &a[j + a_dim1], lda, &c_b12, &a[j + 1 + 
			j * a_dim1], &c__1, (ftnlen)12);
#line 216 "spotf2.f"
		i__2 = *n - j;
#line 216 "spotf2.f"
		d__1 = 1. / ajj;
#line 216 "spotf2.f"
		sscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
#line 217 "spotf2.f"
	    }
#line 218 "spotf2.f"
/* L20: */
#line 218 "spotf2.f"
	}
#line 219 "spotf2.f"
    }
#line 220 "spotf2.f"
    goto L40;

#line 222 "spotf2.f"
L30:
#line 223 "spotf2.f"
    *info = j;

#line 225 "spotf2.f"
L40:
#line 226 "spotf2.f"
    return 0;

/*     End of SPOTF2 */

} /* spotf2_ */

