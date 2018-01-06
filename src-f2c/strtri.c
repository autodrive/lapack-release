#line 1 "strtri.f"
/* strtri.f -- translated by f2c (version 20100827).
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

#line 1 "strtri.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b18 = 1.;
static doublereal c_b22 = -1.;

/* > \brief \b STRTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STRTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strtri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strtri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strtri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRTRI( UPLO, DIAG, N, A, LDA, INFO ) */

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
/* > STRTRI computes the inverse of a real upper or lower triangular */
/* > matrix A. */
/* > */
/* > This is the Level 3 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          = 'N':  A is non-unit triangular; */
/* >          = 'U':  A is unit triangular. */
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
/* >          leading N-by-N upper triangular part of the array A contains */
/* >          the upper triangular matrix, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of the array A contains */
/* >          the lower triangular matrix, and the strictly upper */
/* >          triangular part of A is not referenced.  If DIAG = 'U', the */
/* >          diagonal elements of A are also not referenced and are */
/* >          assumed to be 1. */
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
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular */
/* >               matrix is singular and its inverse can not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int strtri_(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, i__1, i__2[2], i__3, i__4, i__5;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer j, jb, nb, nn;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), strsm_(
	    char *, char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), strti2_(char *, char *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), xerbla_(char 
	    *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical nounit;


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

#line 150 "strtri.f"
    /* Parameter adjustments */
#line 150 "strtri.f"
    a_dim1 = *lda;
#line 150 "strtri.f"
    a_offset = 1 + a_dim1;
#line 150 "strtri.f"
    a -= a_offset;
#line 150 "strtri.f"

#line 150 "strtri.f"
    /* Function Body */
#line 150 "strtri.f"
    *info = 0;
#line 151 "strtri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 152 "strtri.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 153 "strtri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 154 "strtri.f"
	*info = -1;
#line 155 "strtri.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 156 "strtri.f"
	*info = -2;
#line 157 "strtri.f"
    } else if (*n < 0) {
#line 158 "strtri.f"
	*info = -3;
#line 159 "strtri.f"
    } else if (*lda < max(1,*n)) {
#line 160 "strtri.f"
	*info = -5;
#line 161 "strtri.f"
    }
#line 162 "strtri.f"
    if (*info != 0) {
#line 163 "strtri.f"
	i__1 = -(*info);
#line 163 "strtri.f"
	xerbla_("STRTRI", &i__1, (ftnlen)6);
#line 164 "strtri.f"
	return 0;
#line 165 "strtri.f"
    }

/*     Quick return if possible */

#line 169 "strtri.f"
    if (*n == 0) {
#line 169 "strtri.f"
	return 0;
#line 169 "strtri.f"
    }

/*     Check for singularity if non-unit. */

#line 174 "strtri.f"
    if (nounit) {
#line 175 "strtri.f"
	i__1 = *n;
#line 175 "strtri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 176 "strtri.f"
	    if (a[*info + *info * a_dim1] == 0.) {
#line 176 "strtri.f"
		return 0;
#line 176 "strtri.f"
	    }
#line 178 "strtri.f"
/* L10: */
#line 178 "strtri.f"
	}
#line 179 "strtri.f"
	*info = 0;
#line 180 "strtri.f"
    }

/*     Determine the block size for this environment. */

/* Writing concatenation */
#line 184 "strtri.f"
    i__2[0] = 1, a__1[0] = uplo;
#line 184 "strtri.f"
    i__2[1] = 1, a__1[1] = diag;
#line 184 "strtri.f"
    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
#line 184 "strtri.f"
    nb = ilaenv_(&c__1, "STRTRI", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)2);
#line 185 "strtri.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

#line 189 "strtri.f"
	strti2_(uplo, diag, n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);
#line 190 "strtri.f"
    } else {

/*        Use blocked code */

#line 194 "strtri.f"
	if (upper) {

/*           Compute inverse of upper triangular matrix */

#line 198 "strtri.f"
	    i__1 = *n;
#line 198 "strtri.f"
	    i__3 = nb;
#line 198 "strtri.f"
	    for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
/* Computing MIN */
#line 199 "strtri.f"
		i__4 = nb, i__5 = *n - j + 1;
#line 199 "strtri.f"
		jb = min(i__4,i__5);

/*              Compute rows 1:j-1 of current block column */

#line 203 "strtri.f"
		i__4 = j - 1;
#line 203 "strtri.f"
		strmm_("Left", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b18, &a[a_offset], lda, &a[j * a_dim1 + 1], lda, (
			ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 205 "strtri.f"
		i__4 = j - 1;
#line 205 "strtri.f"
		strsm_("Right", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b22, &a[j + j * a_dim1], lda, &a[j * a_dim1 + 1], 
			lda, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)1);

/*              Compute inverse of current diagonal block */

#line 210 "strtri.f"
		strti2_("Upper", diag, &jb, &a[j + j * a_dim1], lda, info, (
			ftnlen)5, (ftnlen)1);
#line 211 "strtri.f"
/* L20: */
#line 211 "strtri.f"
	    }
#line 212 "strtri.f"
	} else {

/*           Compute inverse of lower triangular matrix */

#line 216 "strtri.f"
	    nn = (*n - 1) / nb * nb + 1;
#line 217 "strtri.f"
	    i__3 = -nb;
#line 217 "strtri.f"
	    for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3) {
/* Computing MIN */
#line 218 "strtri.f"
		i__1 = nb, i__4 = *n - j + 1;
#line 218 "strtri.f"
		jb = min(i__1,i__4);
#line 219 "strtri.f"
		if (j + jb <= *n) {

/*                 Compute rows j+jb:n of current block column */

#line 223 "strtri.f"
		    i__1 = *n - j - jb + 1;
#line 223 "strtri.f"
		    strmm_("Left", "Lower", "No transpose", diag, &i__1, &jb, 
			    &c_b18, &a[j + jb + (j + jb) * a_dim1], lda, &a[j 
			    + jb + j * a_dim1], lda, (ftnlen)4, (ftnlen)5, (
			    ftnlen)12, (ftnlen)1);
#line 226 "strtri.f"
		    i__1 = *n - j - jb + 1;
#line 226 "strtri.f"
		    strsm_("Right", "Lower", "No transpose", diag, &i__1, &jb,
			     &c_b22, &a[j + j * a_dim1], lda, &a[j + jb + j * 
			    a_dim1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)12, (
			    ftnlen)1);
#line 229 "strtri.f"
		}

/*              Compute inverse of current diagonal block */

#line 233 "strtri.f"
		strti2_("Lower", diag, &jb, &a[j + j * a_dim1], lda, info, (
			ftnlen)5, (ftnlen)1);
#line 234 "strtri.f"
/* L30: */
#line 234 "strtri.f"
	    }
#line 235 "strtri.f"
	}
#line 236 "strtri.f"
    }

#line 238 "strtri.f"
    return 0;

/*     End of STRTRI */

} /* strtri_ */

