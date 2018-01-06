#line 1 "zpttrf.f"
/* zpttrf.f -- translated by f2c (version 20100827).
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

#line 1 "zpttrf.f"
/* > \brief \b ZPTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPTTRF( N, D, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ) */
/*       COMPLEX*16         E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPTTRF computes the L*D*L**H factorization of a complex Hermitian */
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
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A.  On exit, the n diagonal elements of the diagonal matrix */
/* >          D from the L*D*L**H factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is COMPLEX*16 array, dimension (N-1) */
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

/* > \ingroup complex16PTcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpttrf_(integer *n, doublereal *d__, doublecomplex *e, 
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

#line 128 "zpttrf.f"
    /* Parameter adjustments */
#line 128 "zpttrf.f"
    --e;
#line 128 "zpttrf.f"
    --d__;
#line 128 "zpttrf.f"

#line 128 "zpttrf.f"
    /* Function Body */
#line 128 "zpttrf.f"
    *info = 0;
#line 129 "zpttrf.f"
    if (*n < 0) {
#line 130 "zpttrf.f"
	*info = -1;
#line 131 "zpttrf.f"
	i__1 = -(*info);
#line 131 "zpttrf.f"
	xerbla_("ZPTTRF", &i__1, (ftnlen)6);
#line 132 "zpttrf.f"
	return 0;
#line 133 "zpttrf.f"
    }

/*     Quick return if possible */

#line 137 "zpttrf.f"
    if (*n == 0) {
#line 137 "zpttrf.f"
	return 0;
#line 137 "zpttrf.f"
    }

/*     Compute the L*D*L**H (or U**H *D*U) factorization of A. */

#line 142 "zpttrf.f"
    i4 = (*n - 1) % 4;
#line 143 "zpttrf.f"
    i__1 = i4;
#line 143 "zpttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 144 "zpttrf.f"
	if (d__[i__] <= 0.) {
#line 145 "zpttrf.f"
	    *info = i__;
#line 146 "zpttrf.f"
	    goto L30;
#line 147 "zpttrf.f"
	}
#line 148 "zpttrf.f"
	i__2 = i__;
#line 148 "zpttrf.f"
	eir = e[i__2].r;
#line 149 "zpttrf.f"
	eii = d_imag(&e[i__]);
#line 150 "zpttrf.f"
	f = eir / d__[i__];
#line 151 "zpttrf.f"
	g = eii / d__[i__];
#line 152 "zpttrf.f"
	i__2 = i__;
#line 152 "zpttrf.f"
	z__1.r = f, z__1.i = g;
#line 152 "zpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 153 "zpttrf.f"
	d__[i__ + 1] = d__[i__ + 1] - f * eir - g * eii;
#line 154 "zpttrf.f"
/* L10: */
#line 154 "zpttrf.f"
    }

#line 156 "zpttrf.f"
    i__1 = *n - 4;
#line 156 "zpttrf.f"
    for (i__ = i4 + 1; i__ <= i__1; i__ += 4) {

/*        Drop out of the loop if d(i) <= 0: the matrix is not positive */
/*        definite. */

#line 161 "zpttrf.f"
	if (d__[i__] <= 0.) {
#line 162 "zpttrf.f"
	    *info = i__;
#line 163 "zpttrf.f"
	    goto L30;
#line 164 "zpttrf.f"
	}

/*        Solve for e(i) and d(i+1). */

#line 168 "zpttrf.f"
	i__2 = i__;
#line 168 "zpttrf.f"
	eir = e[i__2].r;
#line 169 "zpttrf.f"
	eii = d_imag(&e[i__]);
#line 170 "zpttrf.f"
	f = eir / d__[i__];
#line 171 "zpttrf.f"
	g = eii / d__[i__];
#line 172 "zpttrf.f"
	i__2 = i__;
#line 172 "zpttrf.f"
	z__1.r = f, z__1.i = g;
#line 172 "zpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 173 "zpttrf.f"
	d__[i__ + 1] = d__[i__ + 1] - f * eir - g * eii;

#line 175 "zpttrf.f"
	if (d__[i__ + 1] <= 0.) {
#line 176 "zpttrf.f"
	    *info = i__ + 1;
#line 177 "zpttrf.f"
	    goto L30;
#line 178 "zpttrf.f"
	}

/*        Solve for e(i+1) and d(i+2). */

#line 182 "zpttrf.f"
	i__2 = i__ + 1;
#line 182 "zpttrf.f"
	eir = e[i__2].r;
#line 183 "zpttrf.f"
	eii = d_imag(&e[i__ + 1]);
#line 184 "zpttrf.f"
	f = eir / d__[i__ + 1];
#line 185 "zpttrf.f"
	g = eii / d__[i__ + 1];
#line 186 "zpttrf.f"
	i__2 = i__ + 1;
#line 186 "zpttrf.f"
	z__1.r = f, z__1.i = g;
#line 186 "zpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 187 "zpttrf.f"
	d__[i__ + 2] = d__[i__ + 2] - f * eir - g * eii;

#line 189 "zpttrf.f"
	if (d__[i__ + 2] <= 0.) {
#line 190 "zpttrf.f"
	    *info = i__ + 2;
#line 191 "zpttrf.f"
	    goto L30;
#line 192 "zpttrf.f"
	}

/*        Solve for e(i+2) and d(i+3). */

#line 196 "zpttrf.f"
	i__2 = i__ + 2;
#line 196 "zpttrf.f"
	eir = e[i__2].r;
#line 197 "zpttrf.f"
	eii = d_imag(&e[i__ + 2]);
#line 198 "zpttrf.f"
	f = eir / d__[i__ + 2];
#line 199 "zpttrf.f"
	g = eii / d__[i__ + 2];
#line 200 "zpttrf.f"
	i__2 = i__ + 2;
#line 200 "zpttrf.f"
	z__1.r = f, z__1.i = g;
#line 200 "zpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 201 "zpttrf.f"
	d__[i__ + 3] = d__[i__ + 3] - f * eir - g * eii;

#line 203 "zpttrf.f"
	if (d__[i__ + 3] <= 0.) {
#line 204 "zpttrf.f"
	    *info = i__ + 3;
#line 205 "zpttrf.f"
	    goto L30;
#line 206 "zpttrf.f"
	}

/*        Solve for e(i+3) and d(i+4). */

#line 210 "zpttrf.f"
	i__2 = i__ + 3;
#line 210 "zpttrf.f"
	eir = e[i__2].r;
#line 211 "zpttrf.f"
	eii = d_imag(&e[i__ + 3]);
#line 212 "zpttrf.f"
	f = eir / d__[i__ + 3];
#line 213 "zpttrf.f"
	g = eii / d__[i__ + 3];
#line 214 "zpttrf.f"
	i__2 = i__ + 3;
#line 214 "zpttrf.f"
	z__1.r = f, z__1.i = g;
#line 214 "zpttrf.f"
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
#line 215 "zpttrf.f"
	d__[i__ + 4] = d__[i__ + 4] - f * eir - g * eii;
#line 216 "zpttrf.f"
/* L20: */
#line 216 "zpttrf.f"
    }

/*     Check d(n) for positive definiteness. */

#line 220 "zpttrf.f"
    if (d__[*n] <= 0.) {
#line 220 "zpttrf.f"
	*info = *n;
#line 220 "zpttrf.f"
    }

#line 223 "zpttrf.f"
L30:
#line 224 "zpttrf.f"
    return 0;

/*     End of ZPTTRF */

} /* zpttrf_ */

