#line 1 "cptcon.f"
/* cptcon.f -- translated by f2c (version 20100827).
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

#line 1 "cptcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CPTCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPTCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPTCON( N, D, E, ANORM, RCOND, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), RWORK( * ) */
/*       COMPLEX            E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPTCON computes the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian positive definite tridiagonal matrix */
/* > using the factorization A = L*D*L**H or A = U**H*D*U computed by */
/* > CPTTRF. */
/* > */
/* > Norm(inv(A)) is computed by a direct method, and the reciprocal of */
/* > the condition number is computed as */
/* >                  RCOND = 1 / (ANORM * norm(inv(A))). */
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
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the diagonal matrix D from the */
/* >          factorization of A, as computed by CPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N-1) */
/* >          The (n-1) off-diagonal elements of the unit bidiagonal factor */
/* >          U or L from the factorization of A, as computed by CPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is REAL */
/* >          The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the */
/* >          1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
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

/* > \date September 2012 */

/* > \ingroup complexPTcomputational */

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
/* Subroutine */ int cptcon_(integer *n, doublereal *d__, doublecomplex *e, 
	doublereal *anorm, doublereal *rcond, doublereal *rwork, integer *
	info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, ix;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal ainvnm;


/*  -- LAPACK computational routine (version 3.4.2) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 160 "cptcon.f"
    /* Parameter adjustments */
#line 160 "cptcon.f"
    --rwork;
#line 160 "cptcon.f"
    --e;
#line 160 "cptcon.f"
    --d__;
#line 160 "cptcon.f"

#line 160 "cptcon.f"
    /* Function Body */
#line 160 "cptcon.f"
    *info = 0;
#line 161 "cptcon.f"
    if (*n < 0) {
#line 162 "cptcon.f"
	*info = -1;
#line 163 "cptcon.f"
    } else if (*anorm < 0.) {
#line 164 "cptcon.f"
	*info = -4;
#line 165 "cptcon.f"
    }
#line 166 "cptcon.f"
    if (*info != 0) {
#line 167 "cptcon.f"
	i__1 = -(*info);
#line 167 "cptcon.f"
	xerbla_("CPTCON", &i__1, (ftnlen)6);
#line 168 "cptcon.f"
	return 0;
#line 169 "cptcon.f"
    }

/*     Quick return if possible */

#line 173 "cptcon.f"
    *rcond = 0.;
#line 174 "cptcon.f"
    if (*n == 0) {
#line 175 "cptcon.f"
	*rcond = 1.;
#line 176 "cptcon.f"
	return 0;
#line 177 "cptcon.f"
    } else if (*anorm == 0.) {
#line 178 "cptcon.f"
	return 0;
#line 179 "cptcon.f"
    }

/*     Check that D(1:N) is positive. */

#line 183 "cptcon.f"
    i__1 = *n;
#line 183 "cptcon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 184 "cptcon.f"
	if (d__[i__] <= 0.) {
#line 184 "cptcon.f"
	    return 0;
#line 184 "cptcon.f"
	}
#line 186 "cptcon.f"
/* L10: */
#line 186 "cptcon.f"
    }

/*     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */

/*        m(i,j) =  abs(A(i,j)), i = j, */
/*        m(i,j) = -abs(A(i,j)), i .ne. j, */

/*     and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H. */

/*     Solve M(L) * x = e. */

#line 197 "cptcon.f"
    rwork[1] = 1.;
#line 198 "cptcon.f"
    i__1 = *n;
#line 198 "cptcon.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 199 "cptcon.f"
	rwork[i__] = rwork[i__ - 1] * z_abs(&e[i__ - 1]) + 1.;
#line 200 "cptcon.f"
/* L20: */
#line 200 "cptcon.f"
    }

/*     Solve D * M(L)**H * x = b. */

#line 204 "cptcon.f"
    rwork[*n] /= d__[*n];
#line 205 "cptcon.f"
    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 206 "cptcon.f"
	rwork[i__] = rwork[i__] / d__[i__] + rwork[i__ + 1] * z_abs(&e[i__]);
#line 207 "cptcon.f"
/* L30: */
#line 207 "cptcon.f"
    }

/*     Compute AINVNM = max(x(i)), 1<=i<=n. */

#line 211 "cptcon.f"
    ix = isamax_(n, &rwork[1], &c__1);
#line 212 "cptcon.f"
    ainvnm = (d__1 = rwork[ix], abs(d__1));

/*     Compute the reciprocal condition number. */

#line 216 "cptcon.f"
    if (ainvnm != 0.) {
#line 216 "cptcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 216 "cptcon.f"
    }

#line 219 "cptcon.f"
    return 0;

/*     End of CPTCON */

} /* cptcon_ */

