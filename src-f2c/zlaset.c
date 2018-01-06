#line 1 "zlaset.f"
/* zlaset.f -- translated by f2c (version 20100827).
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

#line 1 "zlaset.f"
/* > \brief \b ZLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given val
ues. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLASET + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaset.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaset.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaset.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            LDA, M, N */
/*       COMPLEX*16         ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLASET initializes a 2-D array A to BETA on the diagonal and */
/* > ALPHA on the offdiagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies the part of the matrix A to be set. */
/* >          = 'U':      Upper triangular part is set. The lower triangle */
/* >                      is unchanged. */
/* >          = 'L':      Lower triangular part is set. The upper triangle */
/* >                      is unchanged. */
/* >          Otherwise:  All of the matrix A is set. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          On entry, M specifies the number of rows of A. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          On entry, N specifies the number of columns of A. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >          All the offdiagonal array elements are set to ALPHA. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 */
/* >          All the diagonal array elements are set to BETA. */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the m by n matrix A. */
/* >          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j; */
/* >                   A(i,i) = BETA , 1 <= i <= min(m,n) */
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

/* > \date November 2015 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlaset_(char *uplo, integer *m, integer *n, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, integer *
	lda, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 137 "zlaset.f"
    /* Parameter adjustments */
#line 137 "zlaset.f"
    a_dim1 = *lda;
#line 137 "zlaset.f"
    a_offset = 1 + a_dim1;
#line 137 "zlaset.f"
    a -= a_offset;
#line 137 "zlaset.f"

#line 137 "zlaset.f"
    /* Function Body */
#line 137 "zlaset.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Set the diagonal to BETA and the strictly upper triangular */
/*        part of the array to ALPHA. */

#line 142 "zlaset.f"
	i__1 = *n;
#line 142 "zlaset.f"
	for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
#line 143 "zlaset.f"
	    i__3 = j - 1;
#line 143 "zlaset.f"
	    i__2 = min(i__3,*m);
#line 143 "zlaset.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 144 "zlaset.f"
		i__3 = i__ + j * a_dim1;
#line 144 "zlaset.f"
		a[i__3].r = alpha->r, a[i__3].i = alpha->i;
#line 145 "zlaset.f"
/* L10: */
#line 145 "zlaset.f"
	    }
#line 146 "zlaset.f"
/* L20: */
#line 146 "zlaset.f"
	}
#line 147 "zlaset.f"
	i__1 = min(*n,*m);
#line 147 "zlaset.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 148 "zlaset.f"
	    i__2 = i__ + i__ * a_dim1;
#line 148 "zlaset.f"
	    a[i__2].r = beta->r, a[i__2].i = beta->i;
#line 149 "zlaset.f"
/* L30: */
#line 149 "zlaset.f"
	}

#line 151 "zlaset.f"
    } else if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {

/*        Set the diagonal to BETA and the strictly lower triangular */
/*        part of the array to ALPHA. */

#line 156 "zlaset.f"
	i__1 = min(*m,*n);
#line 156 "zlaset.f"
	for (j = 1; j <= i__1; ++j) {
#line 157 "zlaset.f"
	    i__2 = *m;
#line 157 "zlaset.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 158 "zlaset.f"
		i__3 = i__ + j * a_dim1;
#line 158 "zlaset.f"
		a[i__3].r = alpha->r, a[i__3].i = alpha->i;
#line 159 "zlaset.f"
/* L40: */
#line 159 "zlaset.f"
	    }
#line 160 "zlaset.f"
/* L50: */
#line 160 "zlaset.f"
	}
#line 161 "zlaset.f"
	i__1 = min(*n,*m);
#line 161 "zlaset.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 162 "zlaset.f"
	    i__2 = i__ + i__ * a_dim1;
#line 162 "zlaset.f"
	    a[i__2].r = beta->r, a[i__2].i = beta->i;
#line 163 "zlaset.f"
/* L60: */
#line 163 "zlaset.f"
	}

#line 165 "zlaset.f"
    } else {

/*        Set the array to BETA on the diagonal and ALPHA on the */
/*        offdiagonal. */

#line 170 "zlaset.f"
	i__1 = *n;
#line 170 "zlaset.f"
	for (j = 1; j <= i__1; ++j) {
#line 171 "zlaset.f"
	    i__2 = *m;
#line 171 "zlaset.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 172 "zlaset.f"
		i__3 = i__ + j * a_dim1;
#line 172 "zlaset.f"
		a[i__3].r = alpha->r, a[i__3].i = alpha->i;
#line 173 "zlaset.f"
/* L70: */
#line 173 "zlaset.f"
	    }
#line 174 "zlaset.f"
/* L80: */
#line 174 "zlaset.f"
	}
#line 175 "zlaset.f"
	i__1 = min(*m,*n);
#line 175 "zlaset.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 176 "zlaset.f"
	    i__2 = i__ + i__ * a_dim1;
#line 176 "zlaset.f"
	    a[i__2].r = beta->r, a[i__2].i = beta->i;
#line 177 "zlaset.f"
/* L90: */
#line 177 "zlaset.f"
	}
#line 178 "zlaset.f"
    }

#line 180 "zlaset.f"
    return 0;

/*     End of ZLASET */

} /* zlaset_ */

