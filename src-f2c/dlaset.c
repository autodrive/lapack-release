#line 1 "dlaset.f"
/* dlaset.f -- translated by f2c (version 20100827).
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

#line 1 "dlaset.f"
/* > \brief \b DLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given val
ues. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASET + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaset.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaset.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaset.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            LDA, M, N */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASET initializes an m-by-n matrix A to BETA on the diagonal and */
/* > ALPHA on the offdiagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies the part of the matrix A to be set. */
/* >          = 'U':      Upper triangular part is set; the strictly lower */
/* >                      triangular part of A is not changed. */
/* >          = 'L':      Lower triangular part is set; the strictly upper */
/* >                      triangular part of A is not changed. */
/* >          Otherwise:  All of the matrix A is set. */
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
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION */
/* >          The constant to which the offdiagonal elements are to be set. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION */
/* >          The constant to which the diagonal elements are to be set. */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On exit, the leading m-by-n submatrix of A is set as follows: */
/* > */
/* >          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n, */
/* >          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n, */
/* >          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j, */
/* > */
/* >          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlaset_(char *uplo, integer *m, integer *n, doublereal *
	alpha, doublereal *beta, doublereal *a, integer *lda, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

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

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 141 "dlaset.f"
    /* Parameter adjustments */
#line 141 "dlaset.f"
    a_dim1 = *lda;
#line 141 "dlaset.f"
    a_offset = 1 + a_dim1;
#line 141 "dlaset.f"
    a -= a_offset;
#line 141 "dlaset.f"

#line 141 "dlaset.f"
    /* Function Body */
#line 141 "dlaset.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Set the strictly upper triangular or trapezoidal part of the */
/*        array to ALPHA. */

#line 146 "dlaset.f"
	i__1 = *n;
#line 146 "dlaset.f"
	for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 147 "dlaset.f"
	    i__3 = j - 1;
#line 147 "dlaset.f"
	    i__2 = min(i__3,*m);
#line 147 "dlaset.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 148 "dlaset.f"
		a[i__ + j * a_dim1] = *alpha;
#line 149 "dlaset.f"
/* L10: */
#line 149 "dlaset.f"
	    }
#line 150 "dlaset.f"
/* L20: */
#line 150 "dlaset.f"
	}

#line 152 "dlaset.f"
    } else if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {

/*        Set the strictly lower triangular or trapezoidal part of the */
/*        array to ALPHA. */

#line 157 "dlaset.f"
	i__1 = min(*m,*n);
#line 157 "dlaset.f"
	for (j = 1; j <= i__1; ++j) {
#line 158 "dlaset.f"
	    i__2 = *m;
#line 158 "dlaset.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 159 "dlaset.f"
		a[i__ + j * a_dim1] = *alpha;
#line 160 "dlaset.f"
/* L30: */
#line 160 "dlaset.f"
	    }
#line 161 "dlaset.f"
/* L40: */
#line 161 "dlaset.f"
	}

#line 163 "dlaset.f"
    } else {

/*        Set the leading m-by-n submatrix to ALPHA. */

#line 167 "dlaset.f"
	i__1 = *n;
#line 167 "dlaset.f"
	for (j = 1; j <= i__1; ++j) {
#line 168 "dlaset.f"
	    i__2 = *m;
#line 168 "dlaset.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 169 "dlaset.f"
		a[i__ + j * a_dim1] = *alpha;
#line 170 "dlaset.f"
/* L50: */
#line 170 "dlaset.f"
	    }
#line 171 "dlaset.f"
/* L60: */
#line 171 "dlaset.f"
	}
#line 172 "dlaset.f"
    }

/*     Set the first min(M,N) diagonal elements to BETA. */

#line 176 "dlaset.f"
    i__1 = min(*m,*n);
#line 176 "dlaset.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 177 "dlaset.f"
	a[i__ + i__ * a_dim1] = *beta;
#line 178 "dlaset.f"
/* L70: */
#line 178 "dlaset.f"
    }

#line 180 "dlaset.f"
    return 0;

/*     End of DLASET */

} /* dlaset_ */

