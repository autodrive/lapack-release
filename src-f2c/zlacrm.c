#line 1 "zlacrm.f"
/* zlacrm.f -- translated by f2c (version 20100827).
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

#line 1 "zlacrm.f"
/* Table of constant values */

static doublereal c_b6 = 1.;
static doublereal c_b7 = 0.;

/* > \brief \b ZLACRM multiplies a complex matrix by a square real matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLACRM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacrm.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacrm.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacrm.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, LDB, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), C( LDC, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLACRM performs a very simple matrix-matrix multiplication: */
/* >          C := A * B, */
/* > where A is M by N and complex; B is N by N and real; */
/* > C is M by N and complex. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A and of the matrix C. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns and rows of the matrix B and */
/* >          the number of columns of the matrix C. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          A contains the M by N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >=max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
/* >          B contains the N by N matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >=max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC, N) */
/* >          C contains the M by N matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >=max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*M*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlacrm_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *b, integer *ldb, doublecomplex *c__, 
	integer *ldc, doublereal *rwork)
{
    /* System generated locals */
    integer b_dim1, b_offset, a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, l;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible. */

#line 149 "zlacrm.f"
    /* Parameter adjustments */
#line 149 "zlacrm.f"
    a_dim1 = *lda;
#line 149 "zlacrm.f"
    a_offset = 1 + a_dim1;
#line 149 "zlacrm.f"
    a -= a_offset;
#line 149 "zlacrm.f"
    b_dim1 = *ldb;
#line 149 "zlacrm.f"
    b_offset = 1 + b_dim1;
#line 149 "zlacrm.f"
    b -= b_offset;
#line 149 "zlacrm.f"
    c_dim1 = *ldc;
#line 149 "zlacrm.f"
    c_offset = 1 + c_dim1;
#line 149 "zlacrm.f"
    c__ -= c_offset;
#line 149 "zlacrm.f"
    --rwork;
#line 149 "zlacrm.f"

#line 149 "zlacrm.f"
    /* Function Body */
#line 149 "zlacrm.f"
    if (*m == 0 || *n == 0) {
#line 149 "zlacrm.f"
	return 0;
#line 149 "zlacrm.f"
    }

#line 152 "zlacrm.f"
    i__1 = *n;
#line 152 "zlacrm.f"
    for (j = 1; j <= i__1; ++j) {
#line 153 "zlacrm.f"
	i__2 = *m;
#line 153 "zlacrm.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 154 "zlacrm.f"
	    i__3 = i__ + j * a_dim1;
#line 154 "zlacrm.f"
	    rwork[(j - 1) * *m + i__] = a[i__3].r;
#line 155 "zlacrm.f"
/* L10: */
#line 155 "zlacrm.f"
	}
#line 156 "zlacrm.f"
/* L20: */
#line 156 "zlacrm.f"
    }

#line 158 "zlacrm.f"
    l = *m * *n + 1;
#line 159 "zlacrm.f"
    dgemm_("N", "N", m, n, n, &c_b6, &rwork[1], m, &b[b_offset], ldb, &c_b7, &
	    rwork[l], m, (ftnlen)1, (ftnlen)1);
#line 161 "zlacrm.f"
    i__1 = *n;
#line 161 "zlacrm.f"
    for (j = 1; j <= i__1; ++j) {
#line 162 "zlacrm.f"
	i__2 = *m;
#line 162 "zlacrm.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 163 "zlacrm.f"
	    i__3 = i__ + j * c_dim1;
#line 163 "zlacrm.f"
	    i__4 = l + (j - 1) * *m + i__ - 1;
#line 163 "zlacrm.f"
	    c__[i__3].r = rwork[i__4], c__[i__3].i = 0.;
#line 164 "zlacrm.f"
/* L30: */
#line 164 "zlacrm.f"
	}
#line 165 "zlacrm.f"
/* L40: */
#line 165 "zlacrm.f"
    }

#line 167 "zlacrm.f"
    i__1 = *n;
#line 167 "zlacrm.f"
    for (j = 1; j <= i__1; ++j) {
#line 168 "zlacrm.f"
	i__2 = *m;
#line 168 "zlacrm.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 169 "zlacrm.f"
	    rwork[(j - 1) * *m + i__] = d_imag(&a[i__ + j * a_dim1]);
#line 170 "zlacrm.f"
/* L50: */
#line 170 "zlacrm.f"
	}
#line 171 "zlacrm.f"
/* L60: */
#line 171 "zlacrm.f"
    }
#line 172 "zlacrm.f"
    dgemm_("N", "N", m, n, n, &c_b6, &rwork[1], m, &b[b_offset], ldb, &c_b7, &
	    rwork[l], m, (ftnlen)1, (ftnlen)1);
#line 174 "zlacrm.f"
    i__1 = *n;
#line 174 "zlacrm.f"
    for (j = 1; j <= i__1; ++j) {
#line 175 "zlacrm.f"
	i__2 = *m;
#line 175 "zlacrm.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 176 "zlacrm.f"
	    i__3 = i__ + j * c_dim1;
#line 176 "zlacrm.f"
	    i__4 = i__ + j * c_dim1;
#line 176 "zlacrm.f"
	    d__1 = c__[i__4].r;
#line 176 "zlacrm.f"
	    i__5 = l + (j - 1) * *m + i__ - 1;
#line 176 "zlacrm.f"
	    z__1.r = d__1, z__1.i = rwork[i__5];
#line 176 "zlacrm.f"
	    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 178 "zlacrm.f"
/* L70: */
#line 178 "zlacrm.f"
	}
#line 179 "zlacrm.f"
/* L80: */
#line 179 "zlacrm.f"
    }

#line 181 "zlacrm.f"
    return 0;

/*     End of ZLACRM */

} /* zlacrm_ */

