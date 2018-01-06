#line 1 "spttrf.f"
/* spttrf.f -- translated by f2c (version 20100827).
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

#line 1 "spttrf.f"
/* > \brief \b SPTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPTTRF( N, D, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPTTRF computes the L*D*L**T factorization of a real symmetric */
/* > positive definite tridiagonal matrix A.  The factorization may also */
/* > be regarded as having the form A = U**T*D*U. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A.  On exit, the n diagonal elements of the diagonal matrix */
/* >          D from the L*D*L**T factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A.  On exit, the (n-1) subdiagonal elements of the */
/* >          unit bidiagonal factor L from the L*D*L**T factorization of A. */
/* >          E can also be regarded as the superdiagonal of the unit */
/* >          bidiagonal factor U from the U**T*D*U factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* >          > 0: if INFO = k, the leading minor of order k is not */
/* >               positive definite; if k < N, the factorization could not */
/* >               be completed, while if k = N, the factorization was */
/* >               completed, but D(N) <= 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int spttrf_(integer *n, doublereal *d__, doublereal *e, 
	integer *info)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, i4;
    static doublereal ei;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 126 "spttrf.f"
    /* Parameter adjustments */
#line 126 "spttrf.f"
    --e;
#line 126 "spttrf.f"
    --d__;
#line 126 "spttrf.f"

#line 126 "spttrf.f"
    /* Function Body */
#line 126 "spttrf.f"
    *info = 0;
#line 127 "spttrf.f"
    if (*n < 0) {
#line 128 "spttrf.f"
	*info = -1;
#line 129 "spttrf.f"
	i__1 = -(*info);
#line 129 "spttrf.f"
	xerbla_("SPTTRF", &i__1, (ftnlen)6);
#line 130 "spttrf.f"
	return 0;
#line 131 "spttrf.f"
    }

/*     Quick return if possible */

#line 135 "spttrf.f"
    if (*n == 0) {
#line 135 "spttrf.f"
	return 0;
#line 135 "spttrf.f"
    }

/*     Compute the L*D*L**T (or U**T*D*U) factorization of A. */

#line 140 "spttrf.f"
    i4 = (*n - 1) % 4;
#line 141 "spttrf.f"
    i__1 = i4;
#line 141 "spttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 142 "spttrf.f"
	if (d__[i__] <= 0.) {
#line 143 "spttrf.f"
	    *info = i__;
#line 144 "spttrf.f"
	    goto L30;
#line 145 "spttrf.f"
	}
#line 146 "spttrf.f"
	ei = e[i__];
#line 147 "spttrf.f"
	e[i__] = ei / d__[i__];
#line 148 "spttrf.f"
	d__[i__ + 1] -= e[i__] * ei;
#line 149 "spttrf.f"
/* L10: */
#line 149 "spttrf.f"
    }

#line 151 "spttrf.f"
    i__1 = *n - 4;
#line 151 "spttrf.f"
    for (i__ = i4 + 1; i__ <= i__1; i__ += 4) {

/*        Drop out of the loop if d(i) <= 0: the matrix is not positive */
/*        definite. */

#line 156 "spttrf.f"
	if (d__[i__] <= 0.) {
#line 157 "spttrf.f"
	    *info = i__;
#line 158 "spttrf.f"
	    goto L30;
#line 159 "spttrf.f"
	}

/*        Solve for e(i) and d(i+1). */

#line 163 "spttrf.f"
	ei = e[i__];
#line 164 "spttrf.f"
	e[i__] = ei / d__[i__];
#line 165 "spttrf.f"
	d__[i__ + 1] -= e[i__] * ei;

#line 167 "spttrf.f"
	if (d__[i__ + 1] <= 0.) {
#line 168 "spttrf.f"
	    *info = i__ + 1;
#line 169 "spttrf.f"
	    goto L30;
#line 170 "spttrf.f"
	}

/*        Solve for e(i+1) and d(i+2). */

#line 174 "spttrf.f"
	ei = e[i__ + 1];
#line 175 "spttrf.f"
	e[i__ + 1] = ei / d__[i__ + 1];
#line 176 "spttrf.f"
	d__[i__ + 2] -= e[i__ + 1] * ei;

#line 178 "spttrf.f"
	if (d__[i__ + 2] <= 0.) {
#line 179 "spttrf.f"
	    *info = i__ + 2;
#line 180 "spttrf.f"
	    goto L30;
#line 181 "spttrf.f"
	}

/*        Solve for e(i+2) and d(i+3). */

#line 185 "spttrf.f"
	ei = e[i__ + 2];
#line 186 "spttrf.f"
	e[i__ + 2] = ei / d__[i__ + 2];
#line 187 "spttrf.f"
	d__[i__ + 3] -= e[i__ + 2] * ei;

#line 189 "spttrf.f"
	if (d__[i__ + 3] <= 0.) {
#line 190 "spttrf.f"
	    *info = i__ + 3;
#line 191 "spttrf.f"
	    goto L30;
#line 192 "spttrf.f"
	}

/*        Solve for e(i+3) and d(i+4). */

#line 196 "spttrf.f"
	ei = e[i__ + 3];
#line 197 "spttrf.f"
	e[i__ + 3] = ei / d__[i__ + 3];
#line 198 "spttrf.f"
	d__[i__ + 4] -= e[i__ + 3] * ei;
#line 199 "spttrf.f"
/* L20: */
#line 199 "spttrf.f"
    }

/*     Check d(n) for positive definiteness. */

#line 203 "spttrf.f"
    if (d__[*n] <= 0.) {
#line 203 "spttrf.f"
	*info = *n;
#line 203 "spttrf.f"
    }

#line 206 "spttrf.f"
L30:
#line 207 "spttrf.f"
    return 0;

/*     End of SPTTRF */

} /* spttrf_ */

