#line 1 "sptcon.f"
/* sptcon.f -- translated by f2c (version 20100827).
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

#line 1 "sptcon.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SPTCON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPTCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sptcon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sptcon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sptcon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPTCON( N, D, E, ANORM, RCOND, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       REAL               ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPTCON computes the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric positive definite tridiagonal matrix */
/* > using the factorization A = L*D*L**T or A = U**T*D*U computed by */
/* > SPTTRF. */
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
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the diagonal matrix D from the */
/* >          factorization of A, as computed by SPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          The (n-1) off-diagonal elements of the unit bidiagonal factor */
/* >          U or L from the factorization of A,  as computed by SPTTRF. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
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

/* > \ingroup realPTcomputational */

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
/* Subroutine */ int sptcon_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

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

#line 158 "sptcon.f"
    /* Parameter adjustments */
#line 158 "sptcon.f"
    --work;
#line 158 "sptcon.f"
    --e;
#line 158 "sptcon.f"
    --d__;
#line 158 "sptcon.f"

#line 158 "sptcon.f"
    /* Function Body */
#line 158 "sptcon.f"
    *info = 0;
#line 159 "sptcon.f"
    if (*n < 0) {
#line 160 "sptcon.f"
	*info = -1;
#line 161 "sptcon.f"
    } else if (*anorm < 0.) {
#line 162 "sptcon.f"
	*info = -4;
#line 163 "sptcon.f"
    }
#line 164 "sptcon.f"
    if (*info != 0) {
#line 165 "sptcon.f"
	i__1 = -(*info);
#line 165 "sptcon.f"
	xerbla_("SPTCON", &i__1, (ftnlen)6);
#line 166 "sptcon.f"
	return 0;
#line 167 "sptcon.f"
    }

/*     Quick return if possible */

#line 171 "sptcon.f"
    *rcond = 0.;
#line 172 "sptcon.f"
    if (*n == 0) {
#line 173 "sptcon.f"
	*rcond = 1.;
#line 174 "sptcon.f"
	return 0;
#line 175 "sptcon.f"
    } else if (*anorm == 0.) {
#line 176 "sptcon.f"
	return 0;
#line 177 "sptcon.f"
    }

/*     Check that D(1:N) is positive. */

#line 181 "sptcon.f"
    i__1 = *n;
#line 181 "sptcon.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 182 "sptcon.f"
	if (d__[i__] <= 0.) {
#line 182 "sptcon.f"
	    return 0;
#line 182 "sptcon.f"
	}
#line 184 "sptcon.f"
/* L10: */
#line 184 "sptcon.f"
    }

/*     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */

/*        m(i,j) =  abs(A(i,j)), i = j, */
/*        m(i,j) = -abs(A(i,j)), i .ne. j, */

/*     and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T. */

/*     Solve M(L) * x = e. */

#line 195 "sptcon.f"
    work[1] = 1.;
#line 196 "sptcon.f"
    i__1 = *n;
#line 196 "sptcon.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 197 "sptcon.f"
	work[i__] = work[i__ - 1] * (d__1 = e[i__ - 1], abs(d__1)) + 1.;
#line 198 "sptcon.f"
/* L20: */
#line 198 "sptcon.f"
    }

/*     Solve D * M(L)**T * x = b. */

#line 202 "sptcon.f"
    work[*n] /= d__[*n];
#line 203 "sptcon.f"
    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 204 "sptcon.f"
	work[i__] = work[i__] / d__[i__] + work[i__ + 1] * (d__1 = e[i__], 
		abs(d__1));
#line 205 "sptcon.f"
/* L30: */
#line 205 "sptcon.f"
    }

/*     Compute AINVNM = max(x(i)), 1<=i<=n. */

#line 209 "sptcon.f"
    ix = isamax_(n, &work[1], &c__1);
#line 210 "sptcon.f"
    ainvnm = (d__1 = work[ix], abs(d__1));

/*     Compute the reciprocal condition number. */

#line 214 "sptcon.f"
    if (ainvnm != 0.) {
#line 214 "sptcon.f"
	*rcond = 1. / ainvnm / *anorm;
#line 214 "sptcon.f"
    }

#line 217 "sptcon.f"
    return 0;

/*     End of SPTCON */

} /* sptcon_ */

