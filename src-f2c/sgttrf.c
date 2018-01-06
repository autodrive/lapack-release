#line 1 "sgttrf.f"
/* sgttrf.f -- translated by f2c (version 20100827).
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

#line 1 "sgttrf.f"
/* > \brief \b SGTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGTTRF( N, DL, D, DU, DU2, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               D( * ), DL( * ), DU( * ), DU2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGTTRF computes an LU factorization of a real tridiagonal matrix A */
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
/* >          DL is REAL array, dimension (N-1) */
/* >          On entry, DL must contain the (n-1) sub-diagonal elements of */
/* >          A. */
/* > */
/* >          On exit, DL is overwritten by the (n-1) multipliers that */
/* >          define the matrix L from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, D must contain the diagonal elements of A. */
/* > */
/* >          On exit, D is overwritten by the n diagonal elements of the */
/* >          upper triangular matrix U from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* >          DU is REAL array, dimension (N-1) */
/* >          On entry, DU must contain the (n-1) super-diagonal elements */
/* >          of A. */
/* > */
/* >          On exit, DU is overwritten by the (n-1) elements of the first */
/* >          super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] DU2 */
/* > \verbatim */
/* >          DU2 is REAL array, dimension (N-2) */
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

/* > \ingroup realGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgttrf_(integer *n, doublereal *dl, doublereal *d__, 
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

#line 158 "sgttrf.f"
    /* Parameter adjustments */
#line 158 "sgttrf.f"
    --ipiv;
#line 158 "sgttrf.f"
    --du2;
#line 158 "sgttrf.f"
    --du;
#line 158 "sgttrf.f"
    --d__;
#line 158 "sgttrf.f"
    --dl;
#line 158 "sgttrf.f"

#line 158 "sgttrf.f"
    /* Function Body */
#line 158 "sgttrf.f"
    *info = 0;
#line 159 "sgttrf.f"
    if (*n < 0) {
#line 160 "sgttrf.f"
	*info = -1;
#line 161 "sgttrf.f"
	i__1 = -(*info);
#line 161 "sgttrf.f"
	xerbla_("SGTTRF", &i__1, (ftnlen)6);
#line 162 "sgttrf.f"
	return 0;
#line 163 "sgttrf.f"
    }

/*     Quick return if possible */

#line 167 "sgttrf.f"
    if (*n == 0) {
#line 167 "sgttrf.f"
	return 0;
#line 167 "sgttrf.f"
    }

/*     Initialize IPIV(i) = i and DU2(I) = 0 */

#line 172 "sgttrf.f"
    i__1 = *n;
#line 172 "sgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 173 "sgttrf.f"
	ipiv[i__] = i__;
#line 174 "sgttrf.f"
/* L10: */
#line 174 "sgttrf.f"
    }
#line 175 "sgttrf.f"
    i__1 = *n - 2;
#line 175 "sgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 176 "sgttrf.f"
	du2[i__] = 0.;
#line 177 "sgttrf.f"
/* L20: */
#line 177 "sgttrf.f"
    }

#line 179 "sgttrf.f"
    i__1 = *n - 2;
#line 179 "sgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 180 "sgttrf.f"
	if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {

/*           No row interchange required, eliminate DL(I) */

#line 184 "sgttrf.f"
	    if (d__[i__] != 0.) {
#line 185 "sgttrf.f"
		fact = dl[i__] / d__[i__];
#line 186 "sgttrf.f"
		dl[i__] = fact;
#line 187 "sgttrf.f"
		d__[i__ + 1] -= fact * du[i__];
#line 188 "sgttrf.f"
	    }
#line 189 "sgttrf.f"
	} else {

/*           Interchange rows I and I+1, eliminate DL(I) */

#line 193 "sgttrf.f"
	    fact = d__[i__] / dl[i__];
#line 194 "sgttrf.f"
	    d__[i__] = dl[i__];
#line 195 "sgttrf.f"
	    dl[i__] = fact;
#line 196 "sgttrf.f"
	    temp = du[i__];
#line 197 "sgttrf.f"
	    du[i__] = d__[i__ + 1];
#line 198 "sgttrf.f"
	    d__[i__ + 1] = temp - fact * d__[i__ + 1];
#line 199 "sgttrf.f"
	    du2[i__] = du[i__ + 1];
#line 200 "sgttrf.f"
	    du[i__ + 1] = -fact * du[i__ + 1];
#line 201 "sgttrf.f"
	    ipiv[i__] = i__ + 1;
#line 202 "sgttrf.f"
	}
#line 203 "sgttrf.f"
/* L30: */
#line 203 "sgttrf.f"
    }
#line 204 "sgttrf.f"
    if (*n > 1) {
#line 205 "sgttrf.f"
	i__ = *n - 1;
#line 206 "sgttrf.f"
	if ((d__1 = d__[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {
#line 207 "sgttrf.f"
	    if (d__[i__] != 0.) {
#line 208 "sgttrf.f"
		fact = dl[i__] / d__[i__];
#line 209 "sgttrf.f"
		dl[i__] = fact;
#line 210 "sgttrf.f"
		d__[i__ + 1] -= fact * du[i__];
#line 211 "sgttrf.f"
	    }
#line 212 "sgttrf.f"
	} else {
#line 213 "sgttrf.f"
	    fact = d__[i__] / dl[i__];
#line 214 "sgttrf.f"
	    d__[i__] = dl[i__];
#line 215 "sgttrf.f"
	    dl[i__] = fact;
#line 216 "sgttrf.f"
	    temp = du[i__];
#line 217 "sgttrf.f"
	    du[i__] = d__[i__ + 1];
#line 218 "sgttrf.f"
	    d__[i__ + 1] = temp - fact * d__[i__ + 1];
#line 219 "sgttrf.f"
	    ipiv[i__] = i__ + 1;
#line 220 "sgttrf.f"
	}
#line 221 "sgttrf.f"
    }

/*     Check for a zero on the diagonal of U. */

#line 225 "sgttrf.f"
    i__1 = *n;
#line 225 "sgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 226 "sgttrf.f"
	if (d__[i__] == 0.) {
#line 227 "sgttrf.f"
	    *info = i__;
#line 228 "sgttrf.f"
	    goto L50;
#line 229 "sgttrf.f"
	}
#line 230 "sgttrf.f"
/* L40: */
#line 230 "sgttrf.f"
    }
#line 231 "sgttrf.f"
L50:

#line 233 "sgttrf.f"
    return 0;

/*     End of SGTTRF */

} /* sgttrf_ */

