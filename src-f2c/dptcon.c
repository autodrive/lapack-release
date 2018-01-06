#line 1 "dptcon.f"
/* dptcon.f -- translated by f2c (version 20100827).
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

#line 1 "dptcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DPTCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPTCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dptcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dptcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dptcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPTCON( N, D, E, ANORM, RCOND, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTCON computes the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric positive definite tridiagonal matrix */
/* > using the factorization A = L*D*L**T or A = U**T*D*U computed by */
/* > DPTTRF. */
/* > */
/* > Norm(inv(A)) is computed by a direct method, and the reciprocal of */
/* > the condition number is computed as */
/* >              RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the diagonal matrix D from the */
/* >          factorization of A, as computed by DPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) off-diagonal elements of the unit bidiagonal factor */
/* >          U or L from the factorization of A,  as computed by DPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is DOUBLE PRECISION */
/* >          The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the */
/* >          1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup doublePTcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The method used is described in Nicholas J. Higham, "Efficient */
/* >  Algorithms for Computing the Condition Number of a Tridiagonal */
/* >  Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dptcon_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, ix;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;


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

/*     Test the input arguments. */

#line 158 "dptcon.f"
    /* Parameter adjustments */
#line 158 "dptcon.f"
    --work;
#line 158 "dptcon.f"
    --e;
#line 158 "dptcon.f"
    --d__;
#line 158 "dptcon.f"

#line 158 "dptcon.f"
    /* Function Body */
#line 158 "dptcon.f"
    *info = 0;
#line 159 "dptcon.f"
    if (*n < 0) {
#line 160 "dptcon.f"
	*info = -1;
#line 161 "dptcon.f"
    } else if (*anorm < 0.) {
#line 162 "dptcon.f"
	*info = -4;
#line 163 "dptcon.f"
    }
#line 164 "dptcon.f"
    if (*info != 0) {
#line 165 "dptcon.f"
	i__1 = -(*info);
#line 165 "dptcon.f"
	xerbla_("DPTCON", &i__1, (ftnlen)6);
#line 166 "dptcon.f"
	return 0;
#line 167 "dptcon.f"
    }

/*     Quick return if possible */

#line 171 "dptcon.f"
    *rcond = 0.;
#line 172 "dptcon.f"
    if (*n == 0) {
#line 173 "dptcon.f"
	*rcond = 1.;
#line 174 "dptcon.f"
	return 0;
#line 175 "dptcon.f"
    } else if (*anorm == 0.) {
#line 176 "dptcon.f"
	return 0;
#line 177 "dptcon.f"
    }

/*     Check that D(1:N) is positive. */

#line 181 "dptcon.f"
    i__1 = *n;
#line 181 "dptcon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 182 "dptcon.f"
	if (d__[i__] <= 0.) {
#line 182 "dptcon.f"
	    return 0;
#line 182 "dptcon.f"
	}
#line 184 "dptcon.f"
/* L10: */
#line 184 "dptcon.f"
    }

/*     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */

/*        m(i,j) =  abs(A(i,j)), i = j, */
/*        m(i,j) = -abs(A(i,j)), i .ne. j, */

/*     and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T. */

/*     Solve M(L) * x = e. */

#line 195 "dptcon.f"
    work[1] = 1.;
#line 196 "dptcon.f"
    i__1 = *n;
#line 196 "dptcon.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 197 "dptcon.f"
	work[i__] = work[i__ - 1] * (d__1 = e[i__ - 1], abs(d__1)) + 1.;
#line 198 "dptcon.f"
/* L20: */
#line 198 "dptcon.f"
    }

/*     Solve D * M(L)**T * x = b. */

#line 202 "dptcon.f"
    work[*n] /= d__[*n];
#line 203 "dptcon.f"
    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 204 "dptcon.f"
	work[i__] = work[i__] / d__[i__] + work[i__ + 1] * (d__1 = e[i__], 
		abs(d__1));
#line 205 "dptcon.f"
/* L30: */
#line 205 "dptcon.f"
    }

/*     Compute AINVNM = max(x(i)), 1<=i<=n. */

#line 209 "dptcon.f"
    ix = idamax_(n, &work[1], &c__1);
#line 210 "dptcon.f"
    ainvnm = (d__1 = work[ix], abs(d__1));

/*     Compute the reciprocal condition number. */

#line 214 "dptcon.f"
    if (ainvnm != 0.) {
#line 214 "dptcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 214 "dptcon.f"
    }

#line 217 "dptcon.f"
    return 0;

/*     End of DPTCON */

} /* dptcon_ */

