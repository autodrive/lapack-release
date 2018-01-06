#line 1 "dptts2.f"
/* dptts2.f -- translated by f2c (version 20100827).
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

#line 1 "dptts2.f"
/* > \brief \b DPTTS2 solves a tridiagonal system of the form AX=B using the L D LH factorization computed by 
spttrf. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPTTS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dptts2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dptts2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dptts2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPTTS2( N, NRHS, D, E, B, LDB ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTTS2 solves a tridiagonal system of the form */
/* >    A * X = B */
/* > using the L*D*L**T factorization of A computed by DPTTRF.  D is a */
/* > diagonal matrix specified in the vector D, L is a unit bidiagonal */
/* > matrix whose subdiagonal is specified in the vector E, and X and B */
/* > are N by NRHS matrices. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the tridiagonal matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the diagonal matrix D from the */
/* >          L*D*L**T factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of the unit bidiagonal factor */
/* >          L from the L*D*L**T factorization of A.  E can also be regarded */
/* >          as the superdiagonal of the unit bidiagonal factor U from the */
/* >          factorization A = U**T*D*U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side vectors B for the system of */
/* >          linear equations. */
/* >          On exit, the solution vectors, X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doublePTcomputational */

/*  ===================================================================== */
/* Subroutine */ int dptts2_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 129 "dptts2.f"
    /* Parameter adjustments */
#line 129 "dptts2.f"
    --d__;
#line 129 "dptts2.f"
    --e;
#line 129 "dptts2.f"
    b_dim1 = *ldb;
#line 129 "dptts2.f"
    b_offset = 1 + b_dim1;
#line 129 "dptts2.f"
    b -= b_offset;
#line 129 "dptts2.f"

#line 129 "dptts2.f"
    /* Function Body */
#line 129 "dptts2.f"
    if (*n <= 1) {
#line 130 "dptts2.f"
	if (*n == 1) {
#line 130 "dptts2.f"
	    d__1 = 1. / d__[1];
#line 130 "dptts2.f"
	    dscal_(nrhs, &d__1, &b[b_offset], ldb);
#line 130 "dptts2.f"
	}
#line 132 "dptts2.f"
	return 0;
#line 133 "dptts2.f"
    }

/*     Solve A * X = B using the factorization A = L*D*L**T, */
/*     overwriting each right hand side vector with its solution. */

#line 138 "dptts2.f"
    i__1 = *nrhs;
#line 138 "dptts2.f"
    for (j = 1; j <= i__1; ++j) {

/*           Solve L * x = b. */

#line 142 "dptts2.f"
	i__2 = *n;
#line 142 "dptts2.f"
	for (i__ = 2; i__ <= i__2; ++i__) {
#line 143 "dptts2.f"
	    b[i__ + j * b_dim1] -= b[i__ - 1 + j * b_dim1] * e[i__ - 1];
#line 144 "dptts2.f"
/* L10: */
#line 144 "dptts2.f"
	}

/*           Solve D * L**T * x = b. */

#line 148 "dptts2.f"
	b[*n + j * b_dim1] /= d__[*n];
#line 149 "dptts2.f"
	for (i__ = *n - 1; i__ >= 1; --i__) {
#line 150 "dptts2.f"
	    b[i__ + j * b_dim1] = b[i__ + j * b_dim1] / d__[i__] - b[i__ + 1 
		    + j * b_dim1] * e[i__];
#line 151 "dptts2.f"
/* L20: */
#line 151 "dptts2.f"
	}
#line 152 "dptts2.f"
/* L30: */
#line 152 "dptts2.f"
    }

#line 154 "dptts2.f"
    return 0;

/*     End of DPTTS2 */

} /* dptts2_ */

