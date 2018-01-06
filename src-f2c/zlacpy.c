#line 1 "zlacpy.f"
/* zlacpy.f -- translated by f2c (version 20100827).
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

#line 1 "zlacpy.f"
/* > \brief \b ZLACPY copies all or part of one two-dimensional array to another. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLACPY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacpy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacpy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacpy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            LDA, LDB, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLACPY copies all or part of a two-dimensional matrix A to another */
/* > matrix B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies the part of the matrix A to be copied to B. */
/* >          = 'U':      Upper triangular part */
/* >          = 'L':      Lower triangular part */
/* >          Otherwise:  All of the matrix A */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The m by n matrix A.  If UPLO = 'U', only the upper trapezium */
/* >          is accessed; if UPLO = 'L', only the lower trapezium is */
/* >          accessed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          On exit, B = A in the locations specified by UPLO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlacpy_(char *uplo, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 133 "zlacpy.f"
    /* Parameter adjustments */
#line 133 "zlacpy.f"
    a_dim1 = *lda;
#line 133 "zlacpy.f"
    a_offset = 1 + a_dim1;
#line 133 "zlacpy.f"
    a -= a_offset;
#line 133 "zlacpy.f"
    b_dim1 = *ldb;
#line 133 "zlacpy.f"
    b_offset = 1 + b_dim1;
#line 133 "zlacpy.f"
    b -= b_offset;
#line 133 "zlacpy.f"

#line 133 "zlacpy.f"
    /* Function Body */
#line 133 "zlacpy.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 134 "zlacpy.f"
	i__1 = *n;
#line 134 "zlacpy.f"
	for (j = 1; j <= i__1; ++j) {
#line 135 "zlacpy.f"
	    i__2 = min(j,*m);
#line 135 "zlacpy.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 136 "zlacpy.f"
		i__3 = i__ + j * b_dim1;
#line 136 "zlacpy.f"
		i__4 = i__ + j * a_dim1;
#line 136 "zlacpy.f"
		b[i__3].r = a[i__4].r, b[i__3].i = a[i__4].i;
#line 137 "zlacpy.f"
/* L10: */
#line 137 "zlacpy.f"
	    }
#line 138 "zlacpy.f"
/* L20: */
#line 138 "zlacpy.f"
	}

#line 140 "zlacpy.f"
    } else if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 141 "zlacpy.f"
	i__1 = *n;
#line 141 "zlacpy.f"
	for (j = 1; j <= i__1; ++j) {
#line 142 "zlacpy.f"
	    i__2 = *m;
#line 142 "zlacpy.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 143 "zlacpy.f"
		i__3 = i__ + j * b_dim1;
#line 143 "zlacpy.f"
		i__4 = i__ + j * a_dim1;
#line 143 "zlacpy.f"
		b[i__3].r = a[i__4].r, b[i__3].i = a[i__4].i;
#line 144 "zlacpy.f"
/* L30: */
#line 144 "zlacpy.f"
	    }
#line 145 "zlacpy.f"
/* L40: */
#line 145 "zlacpy.f"
	}

#line 147 "zlacpy.f"
    } else {
#line 148 "zlacpy.f"
	i__1 = *n;
#line 148 "zlacpy.f"
	for (j = 1; j <= i__1; ++j) {
#line 149 "zlacpy.f"
	    i__2 = *m;
#line 149 "zlacpy.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 150 "zlacpy.f"
		i__3 = i__ + j * b_dim1;
#line 150 "zlacpy.f"
		i__4 = i__ + j * a_dim1;
#line 150 "zlacpy.f"
		b[i__3].r = a[i__4].r, b[i__3].i = a[i__4].i;
#line 151 "zlacpy.f"
/* L50: */
#line 151 "zlacpy.f"
	    }
#line 152 "zlacpy.f"
/* L60: */
#line 152 "zlacpy.f"
	}
#line 153 "zlacpy.f"
    }

#line 155 "zlacpy.f"
    return 0;

/*     End of ZLACPY */

} /* zlacpy_ */

