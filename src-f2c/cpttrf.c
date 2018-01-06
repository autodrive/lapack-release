#line 1 "cpttrf.f"
/* cpttrf.f -- translated by f2c (version 20100827).
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

#line 1 "cpttrf.f"
/* > \brief \b CPTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPTTRF( N, D, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ) */
/*       COMPLEX            E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPTTRF computes the L*D*L**H factorization of a complex Hermitian */
/* > positive definite tridiagonal matrix A.  The factorization may also */
/* > be regarded as having the form A = U**H *D*U. */
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
/* >          D from the L*D*L**H factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A.  On exit, the (n-1) subdiagonal elements of the */
/* >          unit bidiagonal factor L from the L*D*L**H factorization of A. */
/* >          E can also be regarded as the superdiagonal of the unit */
/* >          bidiagonal factor U from the U**H *D*U factorization of A. */
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

/* > \ingroup complexPTcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpttrf_(integer *n, doublereal *d__, doublecomplex *e, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static doublereal f, g;
    static integer i__, i4;
    static doublereal eii, eir;
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

#line 128 "cpttrf.f"
    /* Parameter adjustments */
#line 128 "cpttrf.f"
    --e;
#line 128 "cpttrf.f"
    --d__;
#line 128 "cpttrf.f"

#line 128 "cpttrf.f"
    /* Function Body */
#line 128 "cpttrf.f"
    *info = 0;
#line 129 "cpttrf.f"
    if (*n < 0) {
#line 130 "cpttrf.f"
	*info = -1;
#line 131 "cpttrf.f"
	i__1 = -(*info);
#line 131 "cpttrf.f"
	xerbla_("CPTTRF", &i__1, (ftnlen)6);
#line 132 "cpttrf.f"
	return 0;
#line 133 "cpttrf.f"
    }

/*     Quick return if possible */

#line 137 "cpttrf.f"
    if (*n == 0) {
#line 137 "cpttrf.f"
	return 0;
#line 137 "cpttrf.f"
    }

/*     Compute the L*D*L**H (or U**H *D*U) factorization of A. */

#line 142 "cpttrf.f"
    i4 = (*n - 1) % 4;
#line 143 "cpttrf.f"
    i__1 = i4;
#line 143 "cpttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 144 "cpttrf.f"
	if (d__[i__] <= 0.) {
#line 145 "cpttrf.f"
	    *info = i__;
#line 146 "cpttrf.f"
	    goto L20;
#line 147 "cpttrf.f"
	}
#line 148 "cpttrf.f"
	i__2 = i__;
#line 148 "cpttrf.f"
	eir = e[i__2].r;
#line 149 "cpttrf.f"
	eii = d_imag(&e[i__]);
#line 150 "cpttrf.f"
	f = eir / d__[i__];
#line 151 "cpttrf.f"
	g = eii / d__[i__];
#line 152 "cpttrf.f"
	i__2 = i__;
#line 152 "cpttrf.f"
	z__1.r = f, z__1.i = g;
#line 152 "cpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 153 "cpttrf.f"
	d__[i__ + 1] = d__[i__ + 1] - f * eir - g * eii;
#line 154 "cpttrf.f"
/* L10: */
#line 154 "cpttrf.f"
    }

#line 156 "cpttrf.f"
    i__1 = *n - 4;
#line 156 "cpttrf.f"
    for (i__ = i4 + 1; i__ <= i__1; i__ += 4) {

/*        Drop out of the loop if d(i) <= 0: the matrix is not positive */
/*        definite. */

#line 161 "cpttrf.f"
	if (d__[i__] <= 0.) {
#line 162 "cpttrf.f"
	    *info = i__;
#line 163 "cpttrf.f"
	    goto L20;
#line 164 "cpttrf.f"
	}

/*        Solve for e(i) and d(i+1). */

#line 168 "cpttrf.f"
	i__2 = i__;
#line 168 "cpttrf.f"
	eir = e[i__2].r;
#line 169 "cpttrf.f"
	eii = d_imag(&e[i__]);
#line 170 "cpttrf.f"
	f = eir / d__[i__];
#line 171 "cpttrf.f"
	g = eii / d__[i__];
#line 172 "cpttrf.f"
	i__2 = i__;
#line 172 "cpttrf.f"
	z__1.r = f, z__1.i = g;
#line 172 "cpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 173 "cpttrf.f"
	d__[i__ + 1] = d__[i__ + 1] - f * eir - g * eii;

#line 175 "cpttrf.f"
	if (d__[i__ + 1] <= 0.) {
#line 176 "cpttrf.f"
	    *info = i__ + 1;
#line 177 "cpttrf.f"
	    goto L20;
#line 178 "cpttrf.f"
	}

/*        Solve for e(i+1) and d(i+2). */

#line 182 "cpttrf.f"
	i__2 = i__ + 1;
#line 182 "cpttrf.f"
	eir = e[i__2].r;
#line 183 "cpttrf.f"
	eii = d_imag(&e[i__ + 1]);
#line 184 "cpttrf.f"
	f = eir / d__[i__ + 1];
#line 185 "cpttrf.f"
	g = eii / d__[i__ + 1];
#line 186 "cpttrf.f"
	i__2 = i__ + 1;
#line 186 "cpttrf.f"
	z__1.r = f, z__1.i = g;
#line 186 "cpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 187 "cpttrf.f"
	d__[i__ + 2] = d__[i__ + 2] - f * eir - g * eii;

#line 189 "cpttrf.f"
	if (d__[i__ + 2] <= 0.) {
#line 190 "cpttrf.f"
	    *info = i__ + 2;
#line 191 "cpttrf.f"
	    goto L20;
#line 192 "cpttrf.f"
	}

/*        Solve for e(i+2) and d(i+3). */

#line 196 "cpttrf.f"
	i__2 = i__ + 2;
#line 196 "cpttrf.f"
	eir = e[i__2].r;
#line 197 "cpttrf.f"
	eii = d_imag(&e[i__ + 2]);
#line 198 "cpttrf.f"
	f = eir / d__[i__ + 2];
#line 199 "cpttrf.f"
	g = eii / d__[i__ + 2];
#line 200 "cpttrf.f"
	i__2 = i__ + 2;
#line 200 "cpttrf.f"
	z__1.r = f, z__1.i = g;
#line 200 "cpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 201 "cpttrf.f"
	d__[i__ + 3] = d__[i__ + 3] - f * eir - g * eii;

#line 203 "cpttrf.f"
	if (d__[i__ + 3] <= 0.) {
#line 204 "cpttrf.f"
	    *info = i__ + 3;
#line 205 "cpttrf.f"
	    goto L20;
#line 206 "cpttrf.f"
	}

/*        Solve for e(i+3) and d(i+4). */

#line 210 "cpttrf.f"
	i__2 = i__ + 3;
#line 210 "cpttrf.f"
	eir = e[i__2].r;
#line 211 "cpttrf.f"
	eii = d_imag(&e[i__ + 3]);
#line 212 "cpttrf.f"
	f = eir / d__[i__ + 3];
#line 213 "cpttrf.f"
	g = eii / d__[i__ + 3];
#line 214 "cpttrf.f"
	i__2 = i__ + 3;
#line 214 "cpttrf.f"
	z__1.r = f, z__1.i = g;
#line 214 "cpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 215 "cpttrf.f"
	d__[i__ + 4] = d__[i__ + 4] - f * eir - g * eii;
#line 216 "cpttrf.f"
/* L110: */
#line 216 "cpttrf.f"
    }

/*     Check d(n) for positive definiteness. */

#line 220 "cpttrf.f"
    if (d__[*n] <= 0.) {
#line 220 "cpttrf.f"
	*info = *n;
#line 220 "cpttrf.f"
    }

#line 223 "cpttrf.f"
L20:
#line 224 "cpttrf.f"
    return 0;

/*     End of CPTTRF */

} /* cpttrf_ */

