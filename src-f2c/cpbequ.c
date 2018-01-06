#line 1 "cpbequ.f"
/* cpbequ.f -- translated by f2c (version 20100827).
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

#line 1 "cpbequ.f"
/* > \brief \b CPBEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPBEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       REAL               AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               S( * ) */
/*       COMPLEX            AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPBEQU computes row and column scalings intended to equilibrate a */
/* > Hermitian positive definite band matrix A and reduce its condition */
/* > number (with respect to the two-norm).  S contains the scale factors, */
/* > S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with */
/* > elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This */
/* > choice of S puts the condition number of B within a factor N of the */
/* > smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangular of A is stored; */
/* >          = 'L':  Lower triangular of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB,N) */
/* >          The upper or lower triangle of the Hermitian band matrix A, */
/* >          stored in the first KD+1 rows of the array.  The j-th column */
/* >          of A is stored in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array A.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL array, dimension (N) */
/* >          If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* >          SCOND is REAL */
/* >          If INFO = 0, S contains the ratio of the smallest S(i) to */
/* >          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too */
/* >          large nor too small, it is not worth scaling by S. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* >          AMAX is REAL */
/* >          Absolute value of largest matrix element.  If AMAX is very */
/* >          close to overflow or very close to underflow, the matrix */
/* >          should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = i, the i-th diagonal element is nonpositive. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpbequ_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal smin;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

/*     Test the input parameters. */

#line 173 "cpbequ.f"
    /* Parameter adjustments */
#line 173 "cpbequ.f"
    ab_dim1 = *ldab;
#line 173 "cpbequ.f"
    ab_offset = 1 + ab_dim1;
#line 173 "cpbequ.f"
    ab -= ab_offset;
#line 173 "cpbequ.f"
    --s;
#line 173 "cpbequ.f"

#line 173 "cpbequ.f"
    /* Function Body */
#line 173 "cpbequ.f"
    *info = 0;
#line 174 "cpbequ.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 175 "cpbequ.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 176 "cpbequ.f"
	*info = -1;
#line 177 "cpbequ.f"
    } else if (*n < 0) {
#line 178 "cpbequ.f"
	*info = -2;
#line 179 "cpbequ.f"
    } else if (*kd < 0) {
#line 180 "cpbequ.f"
	*info = -3;
#line 181 "cpbequ.f"
    } else if (*ldab < *kd + 1) {
#line 182 "cpbequ.f"
	*info = -5;
#line 183 "cpbequ.f"
    }
#line 184 "cpbequ.f"
    if (*info != 0) {
#line 185 "cpbequ.f"
	i__1 = -(*info);
#line 185 "cpbequ.f"
	xerbla_("CPBEQU", &i__1, (ftnlen)6);
#line 186 "cpbequ.f"
	return 0;
#line 187 "cpbequ.f"
    }

/*     Quick return if possible */

#line 191 "cpbequ.f"
    if (*n == 0) {
#line 192 "cpbequ.f"
	*scond = 1.;
#line 193 "cpbequ.f"
	*amax = 0.;
#line 194 "cpbequ.f"
	return 0;
#line 195 "cpbequ.f"
    }

#line 197 "cpbequ.f"
    if (upper) {
#line 198 "cpbequ.f"
	j = *kd + 1;
#line 199 "cpbequ.f"
    } else {
#line 200 "cpbequ.f"
	j = 1;
#line 201 "cpbequ.f"
    }

/*     Initialize SMIN and AMAX. */

#line 205 "cpbequ.f"
    i__1 = j + ab_dim1;
#line 205 "cpbequ.f"
    s[1] = ab[i__1].r;
#line 206 "cpbequ.f"
    smin = s[1];
#line 207 "cpbequ.f"
    *amax = s[1];

/*     Find the minimum and maximum diagonal elements. */

#line 211 "cpbequ.f"
    i__1 = *n;
#line 211 "cpbequ.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 212 "cpbequ.f"
	i__2 = j + i__ * ab_dim1;
#line 212 "cpbequ.f"
	s[i__] = ab[i__2].r;
/* Computing MIN */
#line 213 "cpbequ.f"
	d__1 = smin, d__2 = s[i__];
#line 213 "cpbequ.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 214 "cpbequ.f"
	d__1 = *amax, d__2 = s[i__];
#line 214 "cpbequ.f"
	*amax = max(d__1,d__2);
#line 215 "cpbequ.f"
/* L10: */
#line 215 "cpbequ.f"
    }

#line 217 "cpbequ.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 221 "cpbequ.f"
	i__1 = *n;
#line 221 "cpbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 222 "cpbequ.f"
	    if (s[i__] <= 0.) {
#line 223 "cpbequ.f"
		*info = i__;
#line 224 "cpbequ.f"
		return 0;
#line 225 "cpbequ.f"
	    }
#line 226 "cpbequ.f"
/* L20: */
#line 226 "cpbequ.f"
	}
#line 227 "cpbequ.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 232 "cpbequ.f"
	i__1 = *n;
#line 232 "cpbequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "cpbequ.f"
	    s[i__] = 1. / sqrt(s[i__]);
#line 234 "cpbequ.f"
/* L30: */
#line 234 "cpbequ.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)) */

#line 238 "cpbequ.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 239 "cpbequ.f"
    }
#line 240 "cpbequ.f"
    return 0;

/*     End of CPBEQU */

} /* cpbequ_ */

