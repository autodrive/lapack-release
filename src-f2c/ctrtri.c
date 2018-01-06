#line 1 "ctrtri.f"
/* ctrtri.f -- translated by f2c (version 20100827).
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

#line 1 "ctrtri.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CTRTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTRTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrtri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrtri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrtri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRTRI( UPLO, DIAG, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
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
/* > CTRTRI computes the inverse of a complex upper or lower triangular */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ctrtri_(char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, integer *info, ftnlen uplo_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, i__1, i__2, i__3[2], i__4, i__5;
    doublecomplex z__1;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer j, jb, nb, nn;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    ctrsm_(char *, char *, char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int ctrti2_(char *, char *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
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

#line 151 "ctrtri.f"
    /* Parameter adjustments */
#line 151 "ctrtri.f"
    a_dim1 = *lda;
#line 151 "ctrtri.f"
    a_offset = 1 + a_dim1;
#line 151 "ctrtri.f"
    a -= a_offset;
#line 151 "ctrtri.f"

#line 151 "ctrtri.f"
    /* Function Body */
#line 151 "ctrtri.f"
    *info = 0;
#line 152 "ctrtri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 153 "ctrtri.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 154 "ctrtri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 155 "ctrtri.f"
	*info = -1;
#line 156 "ctrtri.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 157 "ctrtri.f"
	*info = -2;
#line 158 "ctrtri.f"
    } else if (*n < 0) {
#line 159 "ctrtri.f"
	*info = -3;
#line 160 "ctrtri.f"
    } else if (*lda < max(1,*n)) {
#line 161 "ctrtri.f"
	*info = -5;
#line 162 "ctrtri.f"
    }
#line 163 "ctrtri.f"
    if (*info != 0) {
#line 164 "ctrtri.f"
	i__1 = -(*info);
#line 164 "ctrtri.f"
	xerbla_("CTRTRI", &i__1, (ftnlen)6);
#line 165 "ctrtri.f"
	return 0;
#line 166 "ctrtri.f"
    }

/*     Quick return if possible */

#line 170 "ctrtri.f"
    if (*n == 0) {
#line 170 "ctrtri.f"
	return 0;
#line 170 "ctrtri.f"
    }

/*     Check for singularity if non-unit. */

#line 175 "ctrtri.f"
    if (nounit) {
#line 176 "ctrtri.f"
	i__1 = *n;
#line 176 "ctrtri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 177 "ctrtri.f"
	    i__2 = *info + *info * a_dim1;
#line 177 "ctrtri.f"
	    if (a[i__2].r == 0. && a[i__2].i == 0.) {
#line 177 "ctrtri.f"
		return 0;
#line 177 "ctrtri.f"
	    }
#line 179 "ctrtri.f"
/* L10: */
#line 179 "ctrtri.f"
	}
#line 180 "ctrtri.f"
	*info = 0;
#line 181 "ctrtri.f"
    }

/*     Determine the block size for this environment. */

/* Writing concatenation */
#line 185 "ctrtri.f"
    i__3[0] = 1, a__1[0] = uplo;
#line 185 "ctrtri.f"
    i__3[1] = 1, a__1[1] = diag;
#line 185 "ctrtri.f"
    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 185 "ctrtri.f"
    nb = ilaenv_(&c__1, "CTRTRI", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)2);
#line 186 "ctrtri.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

#line 190 "ctrtri.f"
	ctrti2_(uplo, diag, n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);
#line 191 "ctrtri.f"
    } else {

/*        Use blocked code */

#line 195 "ctrtri.f"
	if (upper) {

/*           Compute inverse of upper triangular matrix */

#line 199 "ctrtri.f"
	    i__1 = *n;
#line 199 "ctrtri.f"
	    i__2 = nb;
#line 199 "ctrtri.f"
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 200 "ctrtri.f"
		i__4 = nb, i__5 = *n - j + 1;
#line 200 "ctrtri.f"
		jb = min(i__4,i__5);

/*              Compute rows 1:j-1 of current block column */

#line 204 "ctrtri.f"
		i__4 = j - 1;
#line 204 "ctrtri.f"
		ctrmm_("Left", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b1, &a[a_offset], lda, &a[j * a_dim1 + 1], lda, (
			ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 206 "ctrtri.f"
		i__4 = j - 1;
#line 206 "ctrtri.f"
		z__1.r = -1., z__1.i = -0.;
#line 206 "ctrtri.f"
		ctrsm_("Right", "Upper", "No transpose", diag, &i__4, &jb, &
			z__1, &a[j + j * a_dim1], lda, &a[j * a_dim1 + 1], 
			lda, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)1);

/*              Compute inverse of current diagonal block */

#line 211 "ctrtri.f"
		ctrti2_("Upper", diag, &jb, &a[j + j * a_dim1], lda, info, (
			ftnlen)5, (ftnlen)1);
#line 212 "ctrtri.f"
/* L20: */
#line 212 "ctrtri.f"
	    }
#line 213 "ctrtri.f"
	} else {

/*           Compute inverse of lower triangular matrix */

#line 217 "ctrtri.f"
	    nn = (*n - 1) / nb * nb + 1;
#line 218 "ctrtri.f"
	    i__2 = -nb;
#line 218 "ctrtri.f"
	    for (j = nn; i__2 < 0 ? j >= 1 : j <= 1; j += i__2) {
/* Computing MIN */
#line 219 "ctrtri.f"
		i__1 = nb, i__4 = *n - j + 1;
#line 219 "ctrtri.f"
		jb = min(i__1,i__4);
#line 220 "ctrtri.f"
		if (j + jb <= *n) {

/*                 Compute rows j+jb:n of current block column */

#line 224 "ctrtri.f"
		    i__1 = *n - j - jb + 1;
#line 224 "ctrtri.f"
		    ctrmm_("Left", "Lower", "No transpose", diag, &i__1, &jb, 
			    &c_b1, &a[j + jb + (j + jb) * a_dim1], lda, &a[j 
			    + jb + j * a_dim1], lda, (ftnlen)4, (ftnlen)5, (
			    ftnlen)12, (ftnlen)1);
#line 227 "ctrtri.f"
		    i__1 = *n - j - jb + 1;
#line 227 "ctrtri.f"
		    z__1.r = -1., z__1.i = -0.;
#line 227 "ctrtri.f"
		    ctrsm_("Right", "Lower", "No transpose", diag, &i__1, &jb,
			     &z__1, &a[j + j * a_dim1], lda, &a[j + jb + j * 
			    a_dim1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)12, (
			    ftnlen)1);
#line 230 "ctrtri.f"
		}

/*              Compute inverse of current diagonal block */

#line 234 "ctrtri.f"
		ctrti2_("Lower", diag, &jb, &a[j + j * a_dim1], lda, info, (
			ftnlen)5, (ftnlen)1);
#line 235 "ctrtri.f"
/* L30: */
#line 235 "ctrtri.f"
	    }
#line 236 "ctrtri.f"
	}
#line 237 "ctrtri.f"
    }

#line 239 "ctrtri.f"
    return 0;

/*     End of CTRTRI */

} /* ctrtri_ */

