#line 1 "ztpttr.f"
/* ztpttr.f -- translated by f2c (version 20100827).
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

#line 1 "ztpttr.f"
/* > \brief \b ZTPTTR copies a triangular matrix from the standard packed format (TP) to the standard full for
mat (TR). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTPTTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpttr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpttr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpttr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPTTR( UPLO, N, AP, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N, LDA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPTTR copies a triangular matrix A from standard packed format (TP) */
/* > to standard full format (TR). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular. */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension ( N*(N+1)/2 ), */
/* >          On entry, the upper or lower triangular matrix A, packed */
/* >          columnwise in a linear array. The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, N ) */
/* >          On exit, the triangular matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
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
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ztpttr_(char *uplo, integer *n, doublecomplex *ap, 
	doublecomplex *a, integer *lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 139 "ztpttr.f"
    /* Parameter adjustments */
#line 139 "ztpttr.f"
    --ap;
#line 139 "ztpttr.f"
    a_dim1 = *lda;
#line 139 "ztpttr.f"
    a_offset = 1 + a_dim1;
#line 139 "ztpttr.f"
    a -= a_offset;
#line 139 "ztpttr.f"

#line 139 "ztpttr.f"
    /* Function Body */
#line 139 "ztpttr.f"
    *info = 0;
#line 140 "ztpttr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 141 "ztpttr.f"
    if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 142 "ztpttr.f"
	*info = -1;
#line 143 "ztpttr.f"
    } else if (*n < 0) {
#line 144 "ztpttr.f"
	*info = -2;
#line 145 "ztpttr.f"
    } else if (*lda < max(1,*n)) {
#line 146 "ztpttr.f"
	*info = -5;
#line 147 "ztpttr.f"
    }
#line 148 "ztpttr.f"
    if (*info != 0) {
#line 149 "ztpttr.f"
	i__1 = -(*info);
#line 149 "ztpttr.f"
	xerbla_("ZTPTTR", &i__1, (ftnlen)6);
#line 150 "ztpttr.f"
	return 0;
#line 151 "ztpttr.f"
    }

#line 153 "ztpttr.f"
    if (lower) {
#line 154 "ztpttr.f"
	k = 0;
#line 155 "ztpttr.f"
	i__1 = *n;
#line 155 "ztpttr.f"
	for (j = 1; j <= i__1; ++j) {
#line 156 "ztpttr.f"
	    i__2 = *n;
#line 156 "ztpttr.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 157 "ztpttr.f"
		++k;
#line 158 "ztpttr.f"
		i__3 = i__ + j * a_dim1;
#line 158 "ztpttr.f"
		i__4 = k;
#line 158 "ztpttr.f"
		a[i__3].r = ap[i__4].r, a[i__3].i = ap[i__4].i;
#line 159 "ztpttr.f"
	    }
#line 160 "ztpttr.f"
	}
#line 161 "ztpttr.f"
    } else {
#line 162 "ztpttr.f"
	k = 0;
#line 163 "ztpttr.f"
	i__1 = *n;
#line 163 "ztpttr.f"
	for (j = 1; j <= i__1; ++j) {
#line 164 "ztpttr.f"
	    i__2 = j;
#line 164 "ztpttr.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 165 "ztpttr.f"
		++k;
#line 166 "ztpttr.f"
		i__3 = i__ + j * a_dim1;
#line 166 "ztpttr.f"
		i__4 = k;
#line 166 "ztpttr.f"
		a[i__3].r = ap[i__4].r, a[i__3].i = ap[i__4].i;
#line 167 "ztpttr.f"
	    }
#line 168 "ztpttr.f"
	}
#line 169 "ztpttr.f"
    }


#line 172 "ztpttr.f"
    return 0;

/*     End of ZTPTTR */

} /* ztpttr_ */

