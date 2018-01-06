#line 1 "dpotrf.f"
/* dpotrf.f -- translated by f2c (version 20100827).
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

#line 1 "dpotrf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b13 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b DPOTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPOTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpotrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpotrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpotrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPOTRF computes the Cholesky factorization of a real symmetric */
/* > positive definite matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**T * U,  if UPLO = 'U', or */
/* >    A = L  * L**T,  if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
/* > */
/* > This is the block version of the algorithm, calling Level 3 BLAS. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T. */
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
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the leading minor of order i is not */
/* >                positive definite, and the factorization could not be */
/* >                completed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup doublePOcomputational */

/*  ===================================================================== */
/* Subroutine */ int dpotrf_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer j, jb, nb;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf2_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 148 "dpotrf.f"
    /* Parameter adjustments */
#line 148 "dpotrf.f"
    a_dim1 = *lda;
#line 148 "dpotrf.f"
    a_offset = 1 + a_dim1;
#line 148 "dpotrf.f"
    a -= a_offset;
#line 148 "dpotrf.f"

#line 148 "dpotrf.f"
    /* Function Body */
#line 148 "dpotrf.f"
    *info = 0;
#line 149 "dpotrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 150 "dpotrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 151 "dpotrf.f"
	*info = -1;
#line 152 "dpotrf.f"
    } else if (*n < 0) {
#line 153 "dpotrf.f"
	*info = -2;
#line 154 "dpotrf.f"
    } else if (*lda < max(1,*n)) {
#line 155 "dpotrf.f"
	*info = -4;
#line 156 "dpotrf.f"
    }
#line 157 "dpotrf.f"
    if (*info != 0) {
#line 158 "dpotrf.f"
	i__1 = -(*info);
#line 158 "dpotrf.f"
	xerbla_("DPOTRF", &i__1, (ftnlen)6);
#line 159 "dpotrf.f"
	return 0;
#line 160 "dpotrf.f"
    }

/*     Quick return if possible */

#line 164 "dpotrf.f"
    if (*n == 0) {
#line 164 "dpotrf.f"
	return 0;
#line 164 "dpotrf.f"
    }

/*     Determine the block size for this environment. */

#line 169 "dpotrf.f"
    nb = ilaenv_(&c__1, "DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 170 "dpotrf.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

#line 174 "dpotrf.f"
	dpotrf2_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
#line 175 "dpotrf.f"
    } else {

/*        Use blocked code. */

#line 179 "dpotrf.f"
	if (upper) {

/*           Compute the Cholesky factorization A = U**T*U. */

#line 183 "dpotrf.f"
	    i__1 = *n;
#line 183 "dpotrf.f"
	    i__2 = nb;
#line 183 "dpotrf.f"
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/* Computing MIN */
#line 188 "dpotrf.f"
		i__3 = nb, i__4 = *n - j + 1;
#line 188 "dpotrf.f"
		jb = min(i__3,i__4);
#line 189 "dpotrf.f"
		i__3 = j - 1;
#line 189 "dpotrf.f"
		dsyrk_("Upper", "Transpose", &jb, &i__3, &c_b13, &a[j * 
			a_dim1 + 1], lda, &c_b14, &a[j + j * a_dim1], lda, (
			ftnlen)5, (ftnlen)9);
#line 191 "dpotrf.f"
		dpotrf2_("Upper", &jb, &a[j + j * a_dim1], lda, info, (ftnlen)
			5);
#line 192 "dpotrf.f"
		if (*info != 0) {
#line 192 "dpotrf.f"
		    goto L30;
#line 192 "dpotrf.f"
		}
#line 194 "dpotrf.f"
		if (j + jb <= *n) {

/*                 Compute the current block row. */

#line 198 "dpotrf.f"
		    i__3 = *n - j - jb + 1;
#line 198 "dpotrf.f"
		    i__4 = j - 1;
#line 198 "dpotrf.f"
		    dgemm_("Transpose", "No transpose", &jb, &i__3, &i__4, &
			    c_b13, &a[j * a_dim1 + 1], lda, &a[(j + jb) * 
			    a_dim1 + 1], lda, &c_b14, &a[j + (j + jb) * 
			    a_dim1], lda, (ftnlen)9, (ftnlen)12);
#line 201 "dpotrf.f"
		    i__3 = *n - j - jb + 1;
#line 201 "dpotrf.f"
		    dtrsm_("Left", "Upper", "Transpose", "Non-unit", &jb, &
			    i__3, &c_b14, &a[j + j * a_dim1], lda, &a[j + (j 
			    + jb) * a_dim1], lda, (ftnlen)4, (ftnlen)5, (
			    ftnlen)9, (ftnlen)8);
#line 204 "dpotrf.f"
		}
#line 205 "dpotrf.f"
/* L10: */
#line 205 "dpotrf.f"
	    }

#line 207 "dpotrf.f"
	} else {

/*           Compute the Cholesky factorization A = L*L**T. */

#line 211 "dpotrf.f"
	    i__2 = *n;
#line 211 "dpotrf.f"
	    i__1 = nb;
#line 211 "dpotrf.f"
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/* Computing MIN */
#line 216 "dpotrf.f"
		i__3 = nb, i__4 = *n - j + 1;
#line 216 "dpotrf.f"
		jb = min(i__3,i__4);
#line 217 "dpotrf.f"
		i__3 = j - 1;
#line 217 "dpotrf.f"
		dsyrk_("Lower", "No transpose", &jb, &i__3, &c_b13, &a[j + 
			a_dim1], lda, &c_b14, &a[j + j * a_dim1], lda, (
			ftnlen)5, (ftnlen)12);
#line 219 "dpotrf.f"
		dpotrf2_("Lower", &jb, &a[j + j * a_dim1], lda, info, (ftnlen)
			5);
#line 220 "dpotrf.f"
		if (*info != 0) {
#line 220 "dpotrf.f"
		    goto L30;
#line 220 "dpotrf.f"
		}
#line 222 "dpotrf.f"
		if (j + jb <= *n) {

/*                 Compute the current block column. */

#line 226 "dpotrf.f"
		    i__3 = *n - j - jb + 1;
#line 226 "dpotrf.f"
		    i__4 = j - 1;
#line 226 "dpotrf.f"
		    dgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &
			    c_b13, &a[j + jb + a_dim1], lda, &a[j + a_dim1], 
			    lda, &c_b14, &a[j + jb + j * a_dim1], lda, (
			    ftnlen)12, (ftnlen)9);
#line 229 "dpotrf.f"
		    i__3 = *n - j - jb + 1;
#line 229 "dpotrf.f"
		    dtrsm_("Right", "Lower", "Transpose", "Non-unit", &i__3, &
			    jb, &c_b14, &a[j + j * a_dim1], lda, &a[j + jb + 
			    j * a_dim1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)9,
			     (ftnlen)8);
#line 232 "dpotrf.f"
		}
#line 233 "dpotrf.f"
/* L20: */
#line 233 "dpotrf.f"
	    }
#line 234 "dpotrf.f"
	}
#line 235 "dpotrf.f"
    }
#line 236 "dpotrf.f"
    goto L40;

#line 238 "dpotrf.f"
L30:
#line 239 "dpotrf.f"
    *info = *info + j - 1;

#line 241 "dpotrf.f"
L40:
#line 242 "dpotrf.f"
    return 0;

/*     End of DPOTRF */

} /* dpotrf_ */

