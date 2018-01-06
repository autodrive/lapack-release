#line 1 "zpotrf.f"
/* zpotrf.f -- translated by f2c (version 20100827).
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

#line 1 "zpotrf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b14 = -1.;
static doublereal c_b15 = 1.;

/* > \brief \b ZPOTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPOTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpotrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpotrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpotrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPOTRF computes the Cholesky factorization of a complex Hermitian */
/* > positive definite matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**H * U,  if UPLO = 'U', or */
/* >    A = L  * L**H,  if UPLO = 'L', */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization A = U**H *U or A = L*L**H. */
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

/* > \date November 2011 */

/* > \ingroup complex16POcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpotrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Local variables */
    static integer j, jb, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), zherk_(char *, char *, integer *, 
	    integer *, doublereal *, doublecomplex *, integer *, doublereal *,
	     doublecomplex *, integer *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    zpotf2_(char *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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

#line 149 "zpotrf.f"
    /* Parameter adjustments */
#line 149 "zpotrf.f"
    a_dim1 = *lda;
#line 149 "zpotrf.f"
    a_offset = 1 + a_dim1;
#line 149 "zpotrf.f"
    a -= a_offset;
#line 149 "zpotrf.f"

#line 149 "zpotrf.f"
    /* Function Body */
#line 149 "zpotrf.f"
    *info = 0;
#line 150 "zpotrf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 151 "zpotrf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 152 "zpotrf.f"
	*info = -1;
#line 153 "zpotrf.f"
    } else if (*n < 0) {
#line 154 "zpotrf.f"
	*info = -2;
#line 155 "zpotrf.f"
    } else if (*lda < max(1,*n)) {
#line 156 "zpotrf.f"
	*info = -4;
#line 157 "zpotrf.f"
    }
#line 158 "zpotrf.f"
    if (*info != 0) {
#line 159 "zpotrf.f"
	i__1 = -(*info);
#line 159 "zpotrf.f"
	xerbla_("ZPOTRF", &i__1, (ftnlen)6);
#line 160 "zpotrf.f"
	return 0;
#line 161 "zpotrf.f"
    }

/*     Quick return if possible */

#line 165 "zpotrf.f"
    if (*n == 0) {
#line 165 "zpotrf.f"
	return 0;
#line 165 "zpotrf.f"
    }

/*     Determine the block size for this environment. */

#line 170 "zpotrf.f"
    nb = ilaenv_(&c__1, "ZPOTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 171 "zpotrf.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

#line 175 "zpotrf.f"
	zpotf2_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
#line 176 "zpotrf.f"
    } else {

/*        Use blocked code. */

#line 180 "zpotrf.f"
	if (upper) {

/*           Compute the Cholesky factorization A = U**H *U. */

#line 184 "zpotrf.f"
	    i__1 = *n;
#line 184 "zpotrf.f"
	    i__2 = nb;
#line 184 "zpotrf.f"
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/* Computing MIN */
#line 189 "zpotrf.f"
		i__3 = nb, i__4 = *n - j + 1;
#line 189 "zpotrf.f"
		jb = min(i__3,i__4);
#line 190 "zpotrf.f"
		i__3 = j - 1;
#line 190 "zpotrf.f"
		zherk_("Upper", "Conjugate transpose", &jb, &i__3, &c_b14, &a[
			j * a_dim1 + 1], lda, &c_b15, &a[j + j * a_dim1], lda,
			 (ftnlen)5, (ftnlen)19);
#line 192 "zpotrf.f"
		zpotf2_("Upper", &jb, &a[j + j * a_dim1], lda, info, (ftnlen)
			5);
#line 193 "zpotrf.f"
		if (*info != 0) {
#line 193 "zpotrf.f"
		    goto L30;
#line 193 "zpotrf.f"
		}
#line 195 "zpotrf.f"
		if (j + jb <= *n) {

/*                 Compute the current block row. */

#line 199 "zpotrf.f"
		    i__3 = *n - j - jb + 1;
#line 199 "zpotrf.f"
		    i__4 = j - 1;
#line 199 "zpotrf.f"
		    z__1.r = -1., z__1.i = -0.;
#line 199 "zpotrf.f"
		    zgemm_("Conjugate transpose", "No transpose", &jb, &i__3, 
			    &i__4, &z__1, &a[j * a_dim1 + 1], lda, &a[(j + jb)
			     * a_dim1 + 1], lda, &c_b1, &a[j + (j + jb) * 
			    a_dim1], lda, (ftnlen)19, (ftnlen)12);
#line 203 "zpotrf.f"
		    i__3 = *n - j - jb + 1;
#line 203 "zpotrf.f"
		    ztrsm_("Left", "Upper", "Conjugate transpose", "Non-unit",
			     &jb, &i__3, &c_b1, &a[j + j * a_dim1], lda, &a[j 
			    + (j + jb) * a_dim1], lda, (ftnlen)4, (ftnlen)5, (
			    ftnlen)19, (ftnlen)8);
#line 206 "zpotrf.f"
		}
#line 207 "zpotrf.f"
/* L10: */
#line 207 "zpotrf.f"
	    }

#line 209 "zpotrf.f"
	} else {

/*           Compute the Cholesky factorization A = L*L**H. */

#line 213 "zpotrf.f"
	    i__2 = *n;
#line 213 "zpotrf.f"
	    i__1 = nb;
#line 213 "zpotrf.f"
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/* Computing MIN */
#line 218 "zpotrf.f"
		i__3 = nb, i__4 = *n - j + 1;
#line 218 "zpotrf.f"
		jb = min(i__3,i__4);
#line 219 "zpotrf.f"
		i__3 = j - 1;
#line 219 "zpotrf.f"
		zherk_("Lower", "No transpose", &jb, &i__3, &c_b14, &a[j + 
			a_dim1], lda, &c_b15, &a[j + j * a_dim1], lda, (
			ftnlen)5, (ftnlen)12);
#line 221 "zpotrf.f"
		zpotf2_("Lower", &jb, &a[j + j * a_dim1], lda, info, (ftnlen)
			5);
#line 222 "zpotrf.f"
		if (*info != 0) {
#line 222 "zpotrf.f"
		    goto L30;
#line 222 "zpotrf.f"
		}
#line 224 "zpotrf.f"
		if (j + jb <= *n) {

/*                 Compute the current block column. */

#line 228 "zpotrf.f"
		    i__3 = *n - j - jb + 1;
#line 228 "zpotrf.f"
		    i__4 = j - 1;
#line 228 "zpotrf.f"
		    z__1.r = -1., z__1.i = -0.;
#line 228 "zpotrf.f"
		    zgemm_("No transpose", "Conjugate transpose", &i__3, &jb, 
			    &i__4, &z__1, &a[j + jb + a_dim1], lda, &a[j + 
			    a_dim1], lda, &c_b1, &a[j + jb + j * a_dim1], lda,
			     (ftnlen)12, (ftnlen)19);
#line 232 "zpotrf.f"
		    i__3 = *n - j - jb + 1;
#line 232 "zpotrf.f"
		    ztrsm_("Right", "Lower", "Conjugate transpose", "Non-unit"
			    , &i__3, &jb, &c_b1, &a[j + j * a_dim1], lda, &a[
			    j + jb + j * a_dim1], lda, (ftnlen)5, (ftnlen)5, (
			    ftnlen)19, (ftnlen)8);
#line 235 "zpotrf.f"
		}
#line 236 "zpotrf.f"
/* L20: */
#line 236 "zpotrf.f"
	    }
#line 237 "zpotrf.f"
	}
#line 238 "zpotrf.f"
    }
#line 239 "zpotrf.f"
    goto L40;

#line 241 "zpotrf.f"
L30:
#line 242 "zpotrf.f"
    *info = *info + j - 1;

#line 244 "zpotrf.f"
L40:
#line 245 "zpotrf.f"
    return 0;

/*     End of ZPOTRF */

} /* zpotrf_ */

