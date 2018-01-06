#line 1 "strti2.f"
/* strti2.f -- translated by f2c (version 20100827).
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

#line 1 "strti2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b STRTI2 computes the inverse of a triangular matrix (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STRTI2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strti2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strti2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strti2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRTI2( UPLO, DIAG, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
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
/* > STRTI2 computes the inverse of a real upper or lower triangular */
/* > matrix. */
/* > */
/* > This is the Level 2 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the matrix A is upper or lower triangular. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          Specifies whether or not the matrix A is unit triangular. */
/* >          = 'N':  Non-unit triangular */
/* >          = 'U':  Unit triangular */
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
/* >          On entry, the triangular matrix A.  If UPLO = 'U', the */
/* >          leading n by n upper triangular part of the array A contains */
/* >          the upper triangular matrix, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n by n lower triangular part of the array A contains */
/* >          the lower triangular matrix, and the strictly upper */
/* >          triangular part of A is not referenced.  If DIAG = 'U', the */
/* >          diagonal elements of A are also not referenced and are */
/* >          assumed to be 1. */
/* > */
/* >          On exit, the (triangular) inverse of the original matrix, in */
/* >          the same storage format. */
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
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int strti2_(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer j;
    static doublereal ajj;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int strmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical nounit;


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

#line 151 "strti2.f"
    /* Parameter adjustments */
#line 151 "strti2.f"
    a_dim1 = *lda;
#line 151 "strti2.f"
    a_offset = 1 + a_dim1;
#line 151 "strti2.f"
    a -= a_offset;
#line 151 "strti2.f"

#line 151 "strti2.f"
    /* Function Body */
#line 151 "strti2.f"
    *info = 0;
#line 152 "strti2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 153 "strti2.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 154 "strti2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 155 "strti2.f"
	*info = -1;
#line 156 "strti2.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 157 "strti2.f"
	*info = -2;
#line 158 "strti2.f"
    } else if (*n < 0) {
#line 159 "strti2.f"
	*info = -3;
#line 160 "strti2.f"
    } else if (*lda < max(1,*n)) {
#line 161 "strti2.f"
	*info = -5;
#line 162 "strti2.f"
    }
#line 163 "strti2.f"
    if (*info != 0) {
#line 164 "strti2.f"
	i__1 = -(*info);
#line 164 "strti2.f"
	xerbla_("STRTI2", &i__1, (ftnlen)6);
#line 165 "strti2.f"
	return 0;
#line 166 "strti2.f"
    }

#line 168 "strti2.f"
    if (upper) {

/*        Compute inverse of upper triangular matrix. */

#line 172 "strti2.f"
	i__1 = *n;
#line 172 "strti2.f"
	for (j = 1; j <= i__1; ++j) {
#line 173 "strti2.f"
	    if (nounit) {
#line 174 "strti2.f"
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
#line 175 "strti2.f"
		ajj = -a[j + j * a_dim1];
#line 176 "strti2.f"
	    } else {
#line 177 "strti2.f"
		ajj = -1.;
#line 178 "strti2.f"
	    }

/*           Compute elements 1:j-1 of j-th column. */

#line 182 "strti2.f"
	    i__2 = j - 1;
#line 182 "strti2.f"
	    strmv_("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &
		    a[j * a_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)
		    1);
#line 184 "strti2.f"
	    i__2 = j - 1;
#line 184 "strti2.f"
	    sscal_(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
#line 185 "strti2.f"
/* L10: */
#line 185 "strti2.f"
	}
#line 186 "strti2.f"
    } else {

/*        Compute inverse of lower triangular matrix. */

#line 190 "strti2.f"
	for (j = *n; j >= 1; --j) {
#line 191 "strti2.f"
	    if (nounit) {
#line 192 "strti2.f"
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
#line 193 "strti2.f"
		ajj = -a[j + j * a_dim1];
#line 194 "strti2.f"
	    } else {
#line 195 "strti2.f"
		ajj = -1.;
#line 196 "strti2.f"
	    }
#line 197 "strti2.f"
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

#line 201 "strti2.f"
		i__1 = *n - j;
#line 201 "strti2.f"
		strmv_("Lower", "No transpose", diag, &i__1, &a[j + 1 + (j + 
			1) * a_dim1], lda, &a[j + 1 + j * a_dim1], &c__1, (
			ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 203 "strti2.f"
		i__1 = *n - j;
#line 203 "strti2.f"
		sscal_(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
#line 204 "strti2.f"
	    }
#line 205 "strti2.f"
/* L20: */
#line 205 "strti2.f"
	}
#line 206 "strti2.f"
    }

#line 208 "strti2.f"
    return 0;

/*     End of STRTI2 */

} /* strti2_ */

