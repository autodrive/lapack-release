#line 1 "clarcm.f"
/* clarcm.f -- translated by f2c (version 20100827).
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

#line 1 "clarcm.f"
/* Table of constant values */

static doublereal c_b6 = 1.;
static doublereal c_b7 = 0.;

/* > \brief \b CLARCM copies all or part of a real two-dimensional array to a complex array. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARCM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarcm.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarcm.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarcm.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, LDB, LDC, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), RWORK( * ) */
/*       COMPLEX            B( LDB, * ), C( LDC, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARCM performs a very simple matrix-matrix multiplication: */
/* >          C := A * B, */
/* > where A is M by M and real; B is M by N and complex; */
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
/* >          A is REAL array, dimension (LDA, M) */
/* >          On entry, A contains the M by M matrix A. */
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
/* >          B is COMPLEX array, dimension (LDB, N) */
/* >          On entry, B contains the M by N matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >=max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC, N) */
/* >          On exit, C contains the M by N matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >=max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (2*M*N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clarcm_(integer *m, integer *n, doublereal *a, integer *
	lda, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc,
	 doublereal *rwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, l;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 149 "clarcm.f"
    /* Parameter adjustments */
#line 149 "clarcm.f"
    a_dim1 = *lda;
#line 149 "clarcm.f"
    a_offset = 1 + a_dim1;
#line 149 "clarcm.f"
    a -= a_offset;
#line 149 "clarcm.f"
    b_dim1 = *ldb;
#line 149 "clarcm.f"
    b_offset = 1 + b_dim1;
#line 149 "clarcm.f"
    b -= b_offset;
#line 149 "clarcm.f"
    c_dim1 = *ldc;
#line 149 "clarcm.f"
    c_offset = 1 + c_dim1;
#line 149 "clarcm.f"
    c__ -= c_offset;
#line 149 "clarcm.f"
    --rwork;
#line 149 "clarcm.f"

#line 149 "clarcm.f"
    /* Function Body */
#line 149 "clarcm.f"
    if (*m == 0 || *n == 0) {
#line 149 "clarcm.f"
	return 0;
#line 149 "clarcm.f"
    }

#line 152 "clarcm.f"
    i__1 = *n;
#line 152 "clarcm.f"
    for (j = 1; j <= i__1; ++j) {
#line 153 "clarcm.f"
	i__2 = *m;
#line 153 "clarcm.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 154 "clarcm.f"
	    i__3 = i__ + j * b_dim1;
#line 154 "clarcm.f"
	    rwork[(j - 1) * *m + i__] = b[i__3].r;
#line 155 "clarcm.f"
/* L10: */
#line 155 "clarcm.f"
	}
#line 156 "clarcm.f"
/* L20: */
#line 156 "clarcm.f"
    }

#line 158 "clarcm.f"
    l = *m * *n + 1;
#line 159 "clarcm.f"
    sgemm_("N", "N", m, n, m, &c_b6, &a[a_offset], lda, &rwork[1], m, &c_b7, &
	    rwork[l], m, (ftnlen)1, (ftnlen)1);
#line 161 "clarcm.f"
    i__1 = *n;
#line 161 "clarcm.f"
    for (j = 1; j <= i__1; ++j) {
#line 162 "clarcm.f"
	i__2 = *m;
#line 162 "clarcm.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 163 "clarcm.f"
	    i__3 = i__ + j * c_dim1;
#line 163 "clarcm.f"
	    i__4 = l + (j - 1) * *m + i__ - 1;
#line 163 "clarcm.f"
	    c__[i__3].r = rwork[i__4], c__[i__3].i = 0.;
#line 164 "clarcm.f"
/* L30: */
#line 164 "clarcm.f"
	}
#line 165 "clarcm.f"
/* L40: */
#line 165 "clarcm.f"
    }

#line 167 "clarcm.f"
    i__1 = *n;
#line 167 "clarcm.f"
    for (j = 1; j <= i__1; ++j) {
#line 168 "clarcm.f"
	i__2 = *m;
#line 168 "clarcm.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 169 "clarcm.f"
	    rwork[(j - 1) * *m + i__] = d_imag(&b[i__ + j * b_dim1]);
#line 170 "clarcm.f"
/* L50: */
#line 170 "clarcm.f"
	}
#line 171 "clarcm.f"
/* L60: */
#line 171 "clarcm.f"
    }
#line 172 "clarcm.f"
    sgemm_("N", "N", m, n, m, &c_b6, &a[a_offset], lda, &rwork[1], m, &c_b7, &
	    rwork[l], m, (ftnlen)1, (ftnlen)1);
#line 174 "clarcm.f"
    i__1 = *n;
#line 174 "clarcm.f"
    for (j = 1; j <= i__1; ++j) {
#line 175 "clarcm.f"
	i__2 = *m;
#line 175 "clarcm.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 176 "clarcm.f"
	    i__3 = i__ + j * c_dim1;
#line 176 "clarcm.f"
	    i__4 = i__ + j * c_dim1;
#line 176 "clarcm.f"
	    d__1 = c__[i__4].r;
#line 176 "clarcm.f"
	    i__5 = l + (j - 1) * *m + i__ - 1;
#line 176 "clarcm.f"
	    z__1.r = d__1, z__1.i = rwork[i__5];
#line 176 "clarcm.f"
	    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 178 "clarcm.f"
/* L70: */
#line 178 "clarcm.f"
	}
#line 179 "clarcm.f"
/* L80: */
#line 179 "clarcm.f"
    }

#line 181 "clarcm.f"
    return 0;

/*     End of CLARCM */

} /* clarcm_ */

