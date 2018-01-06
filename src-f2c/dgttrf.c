#line 1 "dgttrf.f"
/* dgttrf.f -- translated by f2c (version 20100827).
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

#line 1 "dgttrf.f"
/* > \brief \b DGTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGTTRF( N, DL, D, DU, DU2, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   D( * ), DL( * ), DU( * ), DU2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGTTRF computes an LU factorization of a real tridiagonal matrix A */
/* > using elimination with partial pivoting and row interchanges. */
/* > */
/* > The factorization has the form */
/* >    A = L * U */
/* > where L is a product of permutation and unit lower bidiagonal */
/* > matrices and U is upper triangular with nonzeros in only the main */
/* > diagonal and first two superdiagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* >          DL is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, DL must contain the (n-1) sub-diagonal elements of */
/* >          A. */
/* > */
/* >          On exit, DL is overwritten by the (n-1) multipliers that */
/* >          define the matrix L from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, D must contain the diagonal elements of A. */
/* > */
/* >          On exit, D is overwritten by the n diagonal elements of the */
/* >          upper triangular matrix U from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* >          DU is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, DU must contain the (n-1) super-diagonal elements */
/* >          of A. */
/* > */
/* >          On exit, DU is overwritten by the (n-1) elements of the first */
/* >          super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] DU2 */
/* > \verbatim */
/* >          DU2 is DOUBLE PRECISION array, dimension (N-2) */
/* >          On exit, DU2 is overwritten by the (n-2) elements of the */
/* >          second super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= n, row i of the matrix was */
/* >          interchanged with row IPIV(i).  IPIV(i) will always be either */
/* >          i or i+1; IPIV(i) = i indicates a row interchange was not */
/* >          required. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -k, the k-th argument had an illegal value */
/* >          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization */
/* >                has been completed, but the factor U is exactly */
/* >                singular, and division by zero will occur if it is used */
/* >                to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgttrf_(integer *n, doublereal *dl, doublereal *d__, 
	doublereal *du, doublereal *du2, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal fact, temp;
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

#line 158 "dgttrf.f"
    /* Parameter adjustments */
#line 158 "dgttrf.f"
    --ipiv;
#line 158 "dgttrf.f"
    --du2;
#line 158 "dgttrf.f"
    --du;
#line 158 "dgttrf.f"
    --d__;
#line 158 "dgttrf.f"
    --dl;
#line 158 "dgttrf.f"

#line 158 "dgttrf.f"
    /* Function Body */
#line 158 "dgttrf.f"
    *info = 0;
#line 159 "dgttrf.f"
    if (*n < 0) {
#line 160 "dgttrf.f"
	*info = -1;
#line 161 "dgttrf.f"
	i__1 = -(*info);
#line 161 "dgttrf.f"
	xerbla_("DGTTRF", &i__1, (ftnlen)6);
#line 162 "dgttrf.f"
	return 0;
#line 163 "dgttrf.f"
    }

/*     Quick return if possible */

#line 167 "dgttrf.f"
    if (*n == 0) {
#line 167 "dgttrf.f"
	return 0;
#line 167 "dgttrf.f"
    }

/*     Initialize IPIV(i) = i and DU2(I) = 0 */

#line 172 "dgttrf.f"
    i__1 = *n;
#line 172 "dgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 173 "dgttrf.f"
	ipiv[i__] = i__;
#line 174 "dgttrf.f"
/* L10: */
#line 174 "dgttrf.f"
    }
#line 175 "dgttrf.f"
    i__1 = *n - 2;
#line 175 "dgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 176 "dgttrf.f"
	du2[i__] = 0.;
#line 177 "dgttrf.f"
/* L20: */
#line 177 "dgttrf.f"
    }

#line 179 "dgttrf.f"
    i__1 = *n - 2;
#line 179 "dgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 180 "dgttrf.f"
	if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {

/*           No row interchange required, eliminate DL(I) */

#line 184 "dgttrf.f"
	    if (d__[i__] != 0.) {
#line 185 "dgttrf.f"
		fact = dl[i__] / d__[i__];
#line 186 "dgttrf.f"
		dl[i__] = fact;
#line 187 "dgttrf.f"
		d__[i__ + 1] -= fact * du[i__];
#line 188 "dgttrf.f"
	    }
#line 189 "dgttrf.f"
	} else {

/*           Interchange rows I and I+1, eliminate DL(I) */

#line 193 "dgttrf.f"
	    fact = d__[i__] / dl[i__];
#line 194 "dgttrf.f"
	    d__[i__] = dl[i__];
#line 195 "dgttrf.f"
	    dl[i__] = fact;
#line 196 "dgttrf.f"
	    temp = du[i__];
#line 197 "dgttrf.f"
	    du[i__] = d__[i__ + 1];
#line 198 "dgttrf.f"
	    d__[i__ + 1] = temp - fact * d__[i__ + 1];
#line 199 "dgttrf.f"
	    du2[i__] = du[i__ + 1];
#line 200 "dgttrf.f"
	    du[i__ + 1] = -fact * du[i__ + 1];
#line 201 "dgttrf.f"
	    ipiv[i__] = i__ + 1;
#line 202 "dgttrf.f"
	}
#line 203 "dgttrf.f"
/* L30: */
#line 203 "dgttrf.f"
    }
#line 204 "dgttrf.f"
    if (*n > 1) {
#line 205 "dgttrf.f"
	i__ = *n - 1;
#line 206 "dgttrf.f"
	if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {
#line 207 "dgttrf.f"
	    if (d__[i__] != 0.) {
#line 208 "dgttrf.f"
		fact = dl[i__] / d__[i__];
#line 209 "dgttrf.f"
		dl[i__] = fact;
#line 210 "dgttrf.f"
		d__[i__ + 1] -= fact * du[i__];
#line 211 "dgttrf.f"
	    }
#line 212 "dgttrf.f"
	} else {
#line 213 "dgttrf.f"
	    fact = d__[i__] / dl[i__];
#line 214 "dgttrf.f"
	    d__[i__] = dl[i__];
#line 215 "dgttrf.f"
	    dl[i__] = fact;
#line 216 "dgttrf.f"
	    temp = du[i__];
#line 217 "dgttrf.f"
	    du[i__] = d__[i__ + 1];
#line 218 "dgttrf.f"
	    d__[i__ + 1] = temp - fact * d__[i__ + 1];
#line 219 "dgttrf.f"
	    ipiv[i__] = i__ + 1;
#line 220 "dgttrf.f"
	}
#line 221 "dgttrf.f"
    }

/*     Check for a zero on the diagonal of U. */

#line 225 "dgttrf.f"
    i__1 = *n;
#line 225 "dgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 226 "dgttrf.f"
	if (d__[i__] == 0.) {
#line 227 "dgttrf.f"
	    *info = i__;
#line 228 "dgttrf.f"
	    goto L50;
#line 229 "dgttrf.f"
	}
#line 230 "dgttrf.f"
/* L40: */
#line 230 "dgttrf.f"
    }
#line 231 "dgttrf.f"
L50:

#line 233 "dgttrf.f"
    return 0;

/*     End of DGTTRF */

} /* dgttrf_ */

