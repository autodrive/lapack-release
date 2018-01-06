#line 1 "dpttrs.f"
/* dpttrs.f -- translated by f2c (version 20100827).
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

#line 1 "dpttrs.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b DPTTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPTTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpttrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpttrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpttrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTTRS solves a tridiagonal system of the form */
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

/* > \date December 2016 */

/* > \ingroup doublePTcomputational */

/*  ===================================================================== */
/* Subroutine */ int dpttrs_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, jb, nb;
    extern /* Subroutine */ int dptts2_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *), xerbla_(char *, integer *,
	     ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 143 "dpttrs.f"
    /* Parameter adjustments */
#line 143 "dpttrs.f"
    --d__;
#line 143 "dpttrs.f"
    --e;
#line 143 "dpttrs.f"
    b_dim1 = *ldb;
#line 143 "dpttrs.f"
    b_offset = 1 + b_dim1;
#line 143 "dpttrs.f"
    b -= b_offset;
#line 143 "dpttrs.f"

#line 143 "dpttrs.f"
    /* Function Body */
#line 143 "dpttrs.f"
    *info = 0;
#line 144 "dpttrs.f"
    if (*n < 0) {
#line 145 "dpttrs.f"
	*info = -1;
#line 146 "dpttrs.f"
    } else if (*nrhs < 0) {
#line 147 "dpttrs.f"
	*info = -2;
#line 148 "dpttrs.f"
    } else if (*ldb < max(1,*n)) {
#line 149 "dpttrs.f"
	*info = -6;
#line 150 "dpttrs.f"
    }
#line 151 "dpttrs.f"
    if (*info != 0) {
#line 152 "dpttrs.f"
	i__1 = -(*info);
#line 152 "dpttrs.f"
	xerbla_("DPTTRS", &i__1, (ftnlen)6);
#line 153 "dpttrs.f"
	return 0;
#line 154 "dpttrs.f"
    }

/*     Quick return if possible */

#line 158 "dpttrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 158 "dpttrs.f"
	return 0;
#line 158 "dpttrs.f"
    }

/*     Determine the number of right-hand sides to solve at a time. */

#line 163 "dpttrs.f"
    if (*nrhs == 1) {
#line 164 "dpttrs.f"
	nb = 1;
#line 165 "dpttrs.f"
    } else {
/* Computing MAX */
#line 166 "dpttrs.f"
	i__1 = 1, i__2 = ilaenv_(&c__1, "DPTTRS", " ", n, nrhs, &c_n1, &c_n1, 
		(ftnlen)6, (ftnlen)1);
#line 166 "dpttrs.f"
	nb = max(i__1,i__2);
#line 167 "dpttrs.f"
    }

#line 169 "dpttrs.f"
    if (nb >= *nrhs) {
#line 170 "dpttrs.f"
	dptts2_(n, nrhs, &d__[1], &e[1], &b[b_offset], ldb);
#line 171 "dpttrs.f"
    } else {
#line 172 "dpttrs.f"
	i__1 = *nrhs;
#line 172 "dpttrs.f"
	i__2 = nb;
#line 172 "dpttrs.f"
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 173 "dpttrs.f"
	    i__3 = *nrhs - j + 1;
#line 173 "dpttrs.f"
	    jb = min(i__3,nb);
#line 174 "dpttrs.f"
	    dptts2_(n, &jb, &d__[1], &e[1], &b[j * b_dim1 + 1], ldb);
#line 175 "dpttrs.f"
/* L10: */
#line 175 "dpttrs.f"
	}
#line 176 "dpttrs.f"
    }

#line 178 "dpttrs.f"
    return 0;

/*     End of DPTTRS */

} /* dpttrs_ */

