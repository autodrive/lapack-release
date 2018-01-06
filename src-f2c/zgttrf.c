#line 1 "zgttrf.f"
/* zgttrf.f -- translated by f2c (version 20100827).
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

#line 1 "zgttrf.f"
/* > \brief \b ZGTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGTTRF( N, DL, D, DU, DU2, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         D( * ), DL( * ), DU( * ), DU2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGTTRF computes an LU factorization of a complex tridiagonal matrix A */
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
/* >          DL is COMPLEX*16 array, dimension (N-1) */
/* >          On entry, DL must contain the (n-1) sub-diagonal elements of */
/* >          A. */
/* > */
/* >          On exit, DL is overwritten by the (n-1) multipliers that */
/* >          define the matrix L from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is COMPLEX*16 array, dimension (N) */
/* >          On entry, D must contain the diagonal elements of A. */
/* > */
/* >          On exit, D is overwritten by the n diagonal elements of the */
/* >          upper triangular matrix U from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* >          DU is COMPLEX*16 array, dimension (N-1) */
/* >          On entry, DU must contain the (n-1) super-diagonal elements */
/* >          of A. */
/* > */
/* >          On exit, DU is overwritten by the (n-1) elements of the first */
/* >          super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] DU2 */
/* > \verbatim */
/* >          DU2 is COMPLEX*16 array, dimension (N-2) */
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

/* > \ingroup complex16GTcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgttrf_(integer *n, doublecomplex *dl, doublecomplex *
	d__, doublecomplex *du, doublecomplex *du2, integer *ipiv, integer *
	info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublecomplex fact, temp;
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 164 "zgttrf.f"
    /* Parameter adjustments */
#line 164 "zgttrf.f"
    --ipiv;
#line 164 "zgttrf.f"
    --du2;
#line 164 "zgttrf.f"
    --du;
#line 164 "zgttrf.f"
    --d__;
#line 164 "zgttrf.f"
    --dl;
#line 164 "zgttrf.f"

#line 164 "zgttrf.f"
    /* Function Body */
#line 164 "zgttrf.f"
    *info = 0;
#line 165 "zgttrf.f"
    if (*n < 0) {
#line 166 "zgttrf.f"
	*info = -1;
#line 167 "zgttrf.f"
	i__1 = -(*info);
#line 167 "zgttrf.f"
	xerbla_("ZGTTRF", &i__1, (ftnlen)6);
#line 168 "zgttrf.f"
	return 0;
#line 169 "zgttrf.f"
    }

/*     Quick return if possible */

#line 173 "zgttrf.f"
    if (*n == 0) {
#line 173 "zgttrf.f"
	return 0;
#line 173 "zgttrf.f"
    }

/*     Initialize IPIV(i) = i and DU2(i) = 0 */

#line 178 "zgttrf.f"
    i__1 = *n;
#line 178 "zgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 179 "zgttrf.f"
	ipiv[i__] = i__;
#line 180 "zgttrf.f"
/* L10: */
#line 180 "zgttrf.f"
    }
#line 181 "zgttrf.f"
    i__1 = *n - 2;
#line 181 "zgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 182 "zgttrf.f"
	i__2 = i__;
#line 182 "zgttrf.f"
	du2[i__2].r = 0., du2[i__2].i = 0.;
#line 183 "zgttrf.f"
/* L20: */
#line 183 "zgttrf.f"
    }

#line 185 "zgttrf.f"
    i__1 = *n - 2;
#line 185 "zgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 186 "zgttrf.f"
	i__2 = i__;
#line 186 "zgttrf.f"
	i__3 = i__;
#line 186 "zgttrf.f"
	if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) >= (d__3 = dl[i__3].r, abs(d__3)) + (d__4 = d_imag(&dl[
		i__]), abs(d__4))) {

/*           No row interchange required, eliminate DL(I) */

#line 190 "zgttrf.f"
	    i__2 = i__;
#line 190 "zgttrf.f"
	    if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), 
		    abs(d__2)) != 0.) {
#line 191 "zgttrf.f"
		z_div(&z__1, &dl[i__], &d__[i__]);
#line 191 "zgttrf.f"
		fact.r = z__1.r, fact.i = z__1.i;
#line 192 "zgttrf.f"
		i__2 = i__;
#line 192 "zgttrf.f"
		dl[i__2].r = fact.r, dl[i__2].i = fact.i;
#line 193 "zgttrf.f"
		i__2 = i__ + 1;
#line 193 "zgttrf.f"
		i__3 = i__ + 1;
#line 193 "zgttrf.f"
		i__4 = i__;
#line 193 "zgttrf.f"
		z__2.r = fact.r * du[i__4].r - fact.i * du[i__4].i, z__2.i = 
			fact.r * du[i__4].i + fact.i * du[i__4].r;
#line 193 "zgttrf.f"
		z__1.r = d__[i__3].r - z__2.r, z__1.i = d__[i__3].i - z__2.i;
#line 193 "zgttrf.f"
		d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
#line 194 "zgttrf.f"
	    }
#line 195 "zgttrf.f"
	} else {

/*           Interchange rows I and I+1, eliminate DL(I) */

#line 199 "zgttrf.f"
	    z_div(&z__1, &d__[i__], &dl[i__]);
#line 199 "zgttrf.f"
	    fact.r = z__1.r, fact.i = z__1.i;
#line 200 "zgttrf.f"
	    i__2 = i__;
#line 200 "zgttrf.f"
	    i__3 = i__;
#line 200 "zgttrf.f"
	    d__[i__2].r = dl[i__3].r, d__[i__2].i = dl[i__3].i;
#line 201 "zgttrf.f"
	    i__2 = i__;
#line 201 "zgttrf.f"
	    dl[i__2].r = fact.r, dl[i__2].i = fact.i;
#line 202 "zgttrf.f"
	    i__2 = i__;
#line 202 "zgttrf.f"
	    temp.r = du[i__2].r, temp.i = du[i__2].i;
#line 203 "zgttrf.f"
	    i__2 = i__;
#line 203 "zgttrf.f"
	    i__3 = i__ + 1;
#line 203 "zgttrf.f"
	    du[i__2].r = d__[i__3].r, du[i__2].i = d__[i__3].i;
#line 204 "zgttrf.f"
	    i__2 = i__ + 1;
#line 204 "zgttrf.f"
	    i__3 = i__ + 1;
#line 204 "zgttrf.f"
	    z__2.r = fact.r * d__[i__3].r - fact.i * d__[i__3].i, z__2.i = 
		    fact.r * d__[i__3].i + fact.i * d__[i__3].r;
#line 204 "zgttrf.f"
	    z__1.r = temp.r - z__2.r, z__1.i = temp.i - z__2.i;
#line 204 "zgttrf.f"
	    d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
#line 205 "zgttrf.f"
	    i__2 = i__;
#line 205 "zgttrf.f"
	    i__3 = i__ + 1;
#line 205 "zgttrf.f"
	    du2[i__2].r = du[i__3].r, du2[i__2].i = du[i__3].i;
#line 206 "zgttrf.f"
	    i__2 = i__ + 1;
#line 206 "zgttrf.f"
	    z__2.r = -fact.r, z__2.i = -fact.i;
#line 206 "zgttrf.f"
	    i__3 = i__ + 1;
#line 206 "zgttrf.f"
	    z__1.r = z__2.r * du[i__3].r - z__2.i * du[i__3].i, z__1.i = 
		    z__2.r * du[i__3].i + z__2.i * du[i__3].r;
#line 206 "zgttrf.f"
	    du[i__2].r = z__1.r, du[i__2].i = z__1.i;
#line 207 "zgttrf.f"
	    ipiv[i__] = i__ + 1;
#line 208 "zgttrf.f"
	}
#line 209 "zgttrf.f"
/* L30: */
#line 209 "zgttrf.f"
    }
#line 210 "zgttrf.f"
    if (*n > 1) {
#line 211 "zgttrf.f"
	i__ = *n - 1;
#line 212 "zgttrf.f"
	i__1 = i__;
#line 212 "zgttrf.f"
	i__2 = i__;
#line 212 "zgttrf.f"
	if ((d__1 = d__[i__1].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) >= (d__3 = dl[i__2].r, abs(d__3)) + (d__4 = d_imag(&dl[
		i__]), abs(d__4))) {
#line 213 "zgttrf.f"
	    i__1 = i__;
#line 213 "zgttrf.f"
	    if ((d__1 = d__[i__1].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), 
		    abs(d__2)) != 0.) {
#line 214 "zgttrf.f"
		z_div(&z__1, &dl[i__], &d__[i__]);
#line 214 "zgttrf.f"
		fact.r = z__1.r, fact.i = z__1.i;
#line 215 "zgttrf.f"
		i__1 = i__;
#line 215 "zgttrf.f"
		dl[i__1].r = fact.r, dl[i__1].i = fact.i;
#line 216 "zgttrf.f"
		i__1 = i__ + 1;
#line 216 "zgttrf.f"
		i__2 = i__ + 1;
#line 216 "zgttrf.f"
		i__3 = i__;
#line 216 "zgttrf.f"
		z__2.r = fact.r * du[i__3].r - fact.i * du[i__3].i, z__2.i = 
			fact.r * du[i__3].i + fact.i * du[i__3].r;
#line 216 "zgttrf.f"
		z__1.r = d__[i__2].r - z__2.r, z__1.i = d__[i__2].i - z__2.i;
#line 216 "zgttrf.f"
		d__[i__1].r = z__1.r, d__[i__1].i = z__1.i;
#line 217 "zgttrf.f"
	    }
#line 218 "zgttrf.f"
	} else {
#line 219 "zgttrf.f"
	    z_div(&z__1, &d__[i__], &dl[i__]);
#line 219 "zgttrf.f"
	    fact.r = z__1.r, fact.i = z__1.i;
#line 220 "zgttrf.f"
	    i__1 = i__;
#line 220 "zgttrf.f"
	    i__2 = i__;
#line 220 "zgttrf.f"
	    d__[i__1].r = dl[i__2].r, d__[i__1].i = dl[i__2].i;
#line 221 "zgttrf.f"
	    i__1 = i__;
#line 221 "zgttrf.f"
	    dl[i__1].r = fact.r, dl[i__1].i = fact.i;
#line 222 "zgttrf.f"
	    i__1 = i__;
#line 222 "zgttrf.f"
	    temp.r = du[i__1].r, temp.i = du[i__1].i;
#line 223 "zgttrf.f"
	    i__1 = i__;
#line 223 "zgttrf.f"
	    i__2 = i__ + 1;
#line 223 "zgttrf.f"
	    du[i__1].r = d__[i__2].r, du[i__1].i = d__[i__2].i;
#line 224 "zgttrf.f"
	    i__1 = i__ + 1;
#line 224 "zgttrf.f"
	    i__2 = i__ + 1;
#line 224 "zgttrf.f"
	    z__2.r = fact.r * d__[i__2].r - fact.i * d__[i__2].i, z__2.i = 
		    fact.r * d__[i__2].i + fact.i * d__[i__2].r;
#line 224 "zgttrf.f"
	    z__1.r = temp.r - z__2.r, z__1.i = temp.i - z__2.i;
#line 224 "zgttrf.f"
	    d__[i__1].r = z__1.r, d__[i__1].i = z__1.i;
#line 225 "zgttrf.f"
	    ipiv[i__] = i__ + 1;
#line 226 "zgttrf.f"
	}
#line 227 "zgttrf.f"
    }

/*     Check for a zero on the diagonal of U. */

#line 231 "zgttrf.f"
    i__1 = *n;
#line 231 "zgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 232 "zgttrf.f"
	i__2 = i__;
#line 232 "zgttrf.f"
	if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) == 0.) {
#line 233 "zgttrf.f"
	    *info = i__;
#line 234 "zgttrf.f"
	    goto L50;
#line 235 "zgttrf.f"
	}
#line 236 "zgttrf.f"
/* L40: */
#line 236 "zgttrf.f"
    }
#line 237 "zgttrf.f"
L50:

#line 239 "zgttrf.f"
    return 0;

/*     End of ZGTTRF */

} /* zgttrf_ */

