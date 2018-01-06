#line 1 "dpttrf.f"
/* dpttrf.f -- translated by f2c (version 20100827).
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

#line 1 "dpttrf.f"
/* > \brief \b DPTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPTTRF( N, D, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTTRF computes the L*D*L**T factorization of a real symmetric */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A.  On exit, the n diagonal elements of the diagonal matrix */
/* >          D from the L*D*L**T factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
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

/* > \date September 2012 */

/* > \ingroup doublePTcomputational */

/*  ===================================================================== */
/* Subroutine */ int dpttrf_(integer *n, doublereal *d__, doublereal *e, 
	integer *info)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, i4;
    static doublereal ei;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 126 "dpttrf.f"
    /* Parameter adjustments */
#line 126 "dpttrf.f"
    --e;
#line 126 "dpttrf.f"
    --d__;
#line 126 "dpttrf.f"

#line 126 "dpttrf.f"
    /* Function Body */
#line 126 "dpttrf.f"
    *info = 0;
#line 127 "dpttrf.f"
    if (*n < 0) {
#line 128 "dpttrf.f"
	*info = -1;
#line 129 "dpttrf.f"
	i__1 = -(*info);
#line 129 "dpttrf.f"
	xerbla_("DPTTRF", &i__1, (ftnlen)6);
#line 130 "dpttrf.f"
	return 0;
#line 131 "dpttrf.f"
    }

/*     Quick return if possible */

#line 135 "dpttrf.f"
    if (*n == 0) {
#line 135 "dpttrf.f"
	return 0;
#line 135 "dpttrf.f"
    }

/*     Compute the L*D*L**T (or U**T*D*U) factorization of A. */

#line 140 "dpttrf.f"
    i4 = (*n - 1) % 4;
#line 141 "dpttrf.f"
    i__1 = i4;
#line 141 "dpttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 142 "dpttrf.f"
	if (d__[i__] <= 0.) {
#line 143 "dpttrf.f"
	    *info = i__;
#line 144 "dpttrf.f"
	    goto L30;
#line 145 "dpttrf.f"
	}
#line 146 "dpttrf.f"
	ei = e[i__];
#line 147 "dpttrf.f"
	e[i__] = ei / d__[i__];
#line 148 "dpttrf.f"
	d__[i__ + 1] -= e[i__] * ei;
#line 149 "dpttrf.f"
/* L10: */
#line 149 "dpttrf.f"
    }

#line 151 "dpttrf.f"
    i__1 = *n - 4;
#line 151 "dpttrf.f"
    for (i__ = i4 + 1; i__ <= i__1; i__ += 4) {

/*        Drop out of the loop if d(i) <= 0: the matrix is not positive */
/*        definite. */

#line 156 "dpttrf.f"
	if (d__[i__] <= 0.) {
#line 157 "dpttrf.f"
	    *info = i__;
#line 158 "dpttrf.f"
	    goto L30;
#line 159 "dpttrf.f"
	}

/*        Solve for e(i) and d(i+1). */

#line 163 "dpttrf.f"
	ei = e[i__];
#line 164 "dpttrf.f"
	e[i__] = ei / d__[i__];
#line 165 "dpttrf.f"
	d__[i__ + 1] -= e[i__] * ei;

#line 167 "dpttrf.f"
	if (d__[i__ + 1] <= 0.) {
#line 168 "dpttrf.f"
	    *info = i__ + 1;
#line 169 "dpttrf.f"
	    goto L30;
#line 170 "dpttrf.f"
	}

/*        Solve for e(i+1) and d(i+2). */

#line 174 "dpttrf.f"
	ei = e[i__ + 1];
#line 175 "dpttrf.f"
	e[i__ + 1] = ei / d__[i__ + 1];
#line 176 "dpttrf.f"
	d__[i__ + 2] -= e[i__ + 1] * ei;

#line 178 "dpttrf.f"
	if (d__[i__ + 2] <= 0.) {
#line 179 "dpttrf.f"
	    *info = i__ + 2;
#line 180 "dpttrf.f"
	    goto L30;
#line 181 "dpttrf.f"
	}

/*        Solve for e(i+2) and d(i+3). */

#line 185 "dpttrf.f"
	ei = e[i__ + 2];
#line 186 "dpttrf.f"
	e[i__ + 2] = ei / d__[i__ + 2];
#line 187 "dpttrf.f"
	d__[i__ + 3] -= e[i__ + 2] * ei;

#line 189 "dpttrf.f"
	if (d__[i__ + 3] <= 0.) {
#line 190 "dpttrf.f"
	    *info = i__ + 3;
#line 191 "dpttrf.f"
	    goto L30;
#line 192 "dpttrf.f"
	}

/*        Solve for e(i+3) and d(i+4). */

#line 196 "dpttrf.f"
	ei = e[i__ + 3];
#line 197 "dpttrf.f"
	e[i__ + 3] = ei / d__[i__ + 3];
#line 198 "dpttrf.f"
	d__[i__ + 4] -= e[i__ + 3] * ei;
#line 199 "dpttrf.f"
/* L20: */
#line 199 "dpttrf.f"
    }

/*     Check d(n) for positive definiteness. */

#line 203 "dpttrf.f"
    if (d__[*n] <= 0.) {
#line 203 "dpttrf.f"
	*info = *n;
#line 203 "dpttrf.f"
    }

#line 206 "dpttrf.f"
L30:
#line 207 "dpttrf.f"
    return 0;

/*     End of DPTTRF */

} /* dpttrf_ */

