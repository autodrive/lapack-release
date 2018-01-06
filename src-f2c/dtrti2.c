#line 1 "dtrti2.f"
/* dtrti2.f -- translated by f2c (version 20100827).
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

#line 1 "dtrti2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DTRTI2 computes the inverse of a triangular matrix (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTRTI2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrti2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrti2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrti2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
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
/* > DTRTI2 computes the inverse of a real upper or lower triangular */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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

/* > \date September 2012 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dtrti2_(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer j;
    static doublereal ajj;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical nounit;


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

#line 151 "dtrti2.f"
    /* Parameter adjustments */
#line 151 "dtrti2.f"
    a_dim1 = *lda;
#line 151 "dtrti2.f"
    a_offset = 1 + a_dim1;
#line 151 "dtrti2.f"
    a -= a_offset;
#line 151 "dtrti2.f"

#line 151 "dtrti2.f"
    /* Function Body */
#line 151 "dtrti2.f"
    *info = 0;
#line 152 "dtrti2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 153 "dtrti2.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 154 "dtrti2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 155 "dtrti2.f"
	*info = -1;
#line 156 "dtrti2.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 157 "dtrti2.f"
	*info = -2;
#line 158 "dtrti2.f"
    } else if (*n < 0) {
#line 159 "dtrti2.f"
	*info = -3;
#line 160 "dtrti2.f"
    } else if (*lda < max(1,*n)) {
#line 161 "dtrti2.f"
	*info = -5;
#line 162 "dtrti2.f"
    }
#line 163 "dtrti2.f"
    if (*info != 0) {
#line 164 "dtrti2.f"
	i__1 = -(*info);
#line 164 "dtrti2.f"
	xerbla_("DTRTI2", &i__1, (ftnlen)6);
#line 165 "dtrti2.f"
	return 0;
#line 166 "dtrti2.f"
    }

#line 168 "dtrti2.f"
    if (upper) {

/*        Compute inverse of upper triangular matrix. */

#line 172 "dtrti2.f"
	i__1 = *n;
#line 172 "dtrti2.f"
	for (j = 1; j <= i__1; ++j) {
#line 173 "dtrti2.f"
	    if (nounit) {
#line 174 "dtrti2.f"
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
#line 175 "dtrti2.f"
		ajj = -a[j + j * a_dim1];
#line 176 "dtrti2.f"
	    } else {
#line 177 "dtrti2.f"
		ajj = -1.;
#line 178 "dtrti2.f"
	    }

/*           Compute elements 1:j-1 of j-th column. */

#line 182 "dtrti2.f"
	    i__2 = j - 1;
#line 182 "dtrti2.f"
	    dtrmv_("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &
		    a[j * a_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)
		    1);
#line 184 "dtrti2.f"
	    i__2 = j - 1;
#line 184 "dtrti2.f"
	    dscal_(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
#line 185 "dtrti2.f"
/* L10: */
#line 185 "dtrti2.f"
	}
#line 186 "dtrti2.f"
    } else {

/*        Compute inverse of lower triangular matrix. */

#line 190 "dtrti2.f"
	for (j = *n; j >= 1; --j) {
#line 191 "dtrti2.f"
	    if (nounit) {
#line 192 "dtrti2.f"
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
#line 193 "dtrti2.f"
		ajj = -a[j + j * a_dim1];
#line 194 "dtrti2.f"
	    } else {
#line 195 "dtrti2.f"
		ajj = -1.;
#line 196 "dtrti2.f"
	    }
#line 197 "dtrti2.f"
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

#line 201 "dtrti2.f"
		i__1 = *n - j;
#line 201 "dtrti2.f"
		dtrmv_("Lower", "No transpose", diag, &i__1, &a[j + 1 + (j + 
			1) * a_dim1], lda, &a[j + 1 + j * a_dim1], &c__1, (
			ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 203 "dtrti2.f"
		i__1 = *n - j;
#line 203 "dtrti2.f"
		dscal_(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
#line 204 "dtrti2.f"
	    }
#line 205 "dtrti2.f"
/* L20: */
#line 205 "dtrti2.f"
	}
#line 206 "dtrti2.f"
    }

#line 208 "dtrti2.f"
    return 0;

/*     End of DTRTI2 */

} /* dtrti2_ */

