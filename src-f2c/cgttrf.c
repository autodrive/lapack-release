#line 1 "cgttrf.f"
/* cgttrf.f -- translated by f2c (version 20100827).
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

#line 1 "cgttrf.f"
/* > \brief \b CGTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGTTRF( N, DL, D, DU, DU2, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            D( * ), DL( * ), DU( * ), DU2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGTTRF computes an LU factorization of a complex tridiagonal matrix A */
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
/* >          DL is COMPLEX array, dimension (N-1) */
/* >          On entry, DL must contain the (n-1) sub-diagonal elements of */
/* >          A. */
/* > */
/* >          On exit, DL is overwritten by the (n-1) multipliers that */
/* >          define the matrix L from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          On entry, D must contain the diagonal elements of A. */
/* > */
/* >          On exit, D is overwritten by the n diagonal elements of the */
/* >          upper triangular matrix U from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* >          DU is COMPLEX array, dimension (N-1) */
/* >          On entry, DU must contain the (n-1) super-diagonal elements */
/* >          of A. */
/* > */
/* >          On exit, DU is overwritten by the (n-1) elements of the first */
/* >          super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] DU2 */
/* > \verbatim */
/* >          DU2 is COMPLEX array, dimension (N-2) */
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

/* > \date December 2016 */

/* > \ingroup complexGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgttrf_(integer *n, doublecomplex *dl, doublecomplex *
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 164 "cgttrf.f"
    /* Parameter adjustments */
#line 164 "cgttrf.f"
    --ipiv;
#line 164 "cgttrf.f"
    --du2;
#line 164 "cgttrf.f"
    --du;
#line 164 "cgttrf.f"
    --d__;
#line 164 "cgttrf.f"
    --dl;
#line 164 "cgttrf.f"

#line 164 "cgttrf.f"
    /* Function Body */
#line 164 "cgttrf.f"
    *info = 0;
#line 165 "cgttrf.f"
    if (*n < 0) {
#line 166 "cgttrf.f"
	*info = -1;
#line 167 "cgttrf.f"
	i__1 = -(*info);
#line 167 "cgttrf.f"
	xerbla_("CGTTRF", &i__1, (ftnlen)6);
#line 168 "cgttrf.f"
	return 0;
#line 169 "cgttrf.f"
    }

/*     Quick return if possible */

#line 173 "cgttrf.f"
    if (*n == 0) {
#line 173 "cgttrf.f"
	return 0;
#line 173 "cgttrf.f"
    }

/*     Initialize IPIV(i) = i and DU2(i) = 0 */

#line 178 "cgttrf.f"
    i__1 = *n;
#line 178 "cgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 179 "cgttrf.f"
	ipiv[i__] = i__;
#line 180 "cgttrf.f"
/* L10: */
#line 180 "cgttrf.f"
    }
#line 181 "cgttrf.f"
    i__1 = *n - 2;
#line 181 "cgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 182 "cgttrf.f"
	i__2 = i__;
#line 182 "cgttrf.f"
	du2[i__2].r = 0., du2[i__2].i = 0.;
#line 183 "cgttrf.f"
/* L20: */
#line 183 "cgttrf.f"
    }

#line 185 "cgttrf.f"
    i__1 = *n - 2;
#line 185 "cgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 186 "cgttrf.f"
	i__2 = i__;
#line 186 "cgttrf.f"
	i__3 = i__;
#line 186 "cgttrf.f"
	if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) >= (d__3 = dl[i__3].r, abs(d__3)) + (d__4 = d_imag(&dl[
		i__]), abs(d__4))) {

/*           No row interchange required, eliminate DL(I) */

#line 190 "cgttrf.f"
	    i__2 = i__;
#line 190 "cgttrf.f"
	    if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), 
		    abs(d__2)) != 0.) {
#line 191 "cgttrf.f"
		z_div(&z__1, &dl[i__], &d__[i__]);
#line 191 "cgttrf.f"
		fact.r = z__1.r, fact.i = z__1.i;
#line 192 "cgttrf.f"
		i__2 = i__;
#line 192 "cgttrf.f"
		dl[i__2].r = fact.r, dl[i__2].i = fact.i;
#line 193 "cgttrf.f"
		i__2 = i__ + 1;
#line 193 "cgttrf.f"
		i__3 = i__ + 1;
#line 193 "cgttrf.f"
		i__4 = i__;
#line 193 "cgttrf.f"
		z__2.r = fact.r * du[i__4].r - fact.i * du[i__4].i, z__2.i = 
			fact.r * du[i__4].i + fact.i * du[i__4].r;
#line 193 "cgttrf.f"
		z__1.r = d__[i__3].r - z__2.r, z__1.i = d__[i__3].i - z__2.i;
#line 193 "cgttrf.f"
		d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
#line 194 "cgttrf.f"
	    }
#line 195 "cgttrf.f"
	} else {

/*           Interchange rows I and I+1, eliminate DL(I) */

#line 199 "cgttrf.f"
	    z_div(&z__1, &d__[i__], &dl[i__]);
#line 199 "cgttrf.f"
	    fact.r = z__1.r, fact.i = z__1.i;
#line 200 "cgttrf.f"
	    i__2 = i__;
#line 200 "cgttrf.f"
	    i__3 = i__;
#line 200 "cgttrf.f"
	    d__[i__2].r = dl[i__3].r, d__[i__2].i = dl[i__3].i;
#line 201 "cgttrf.f"
	    i__2 = i__;
#line 201 "cgttrf.f"
	    dl[i__2].r = fact.r, dl[i__2].i = fact.i;
#line 202 "cgttrf.f"
	    i__2 = i__;
#line 202 "cgttrf.f"
	    temp.r = du[i__2].r, temp.i = du[i__2].i;
#line 203 "cgttrf.f"
	    i__2 = i__;
#line 203 "cgttrf.f"
	    i__3 = i__ + 1;
#line 203 "cgttrf.f"
	    du[i__2].r = d__[i__3].r, du[i__2].i = d__[i__3].i;
#line 204 "cgttrf.f"
	    i__2 = i__ + 1;
#line 204 "cgttrf.f"
	    i__3 = i__ + 1;
#line 204 "cgttrf.f"
	    z__2.r = fact.r * d__[i__3].r - fact.i * d__[i__3].i, z__2.i = 
		    fact.r * d__[i__3].i + fact.i * d__[i__3].r;
#line 204 "cgttrf.f"
	    z__1.r = temp.r - z__2.r, z__1.i = temp.i - z__2.i;
#line 204 "cgttrf.f"
	    d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
#line 205 "cgttrf.f"
	    i__2 = i__;
#line 205 "cgttrf.f"
	    i__3 = i__ + 1;
#line 205 "cgttrf.f"
	    du2[i__2].r = du[i__3].r, du2[i__2].i = du[i__3].i;
#line 206 "cgttrf.f"
	    i__2 = i__ + 1;
#line 206 "cgttrf.f"
	    z__2.r = -fact.r, z__2.i = -fact.i;
#line 206 "cgttrf.f"
	    i__3 = i__ + 1;
#line 206 "cgttrf.f"
	    z__1.r = z__2.r * du[i__3].r - z__2.i * du[i__3].i, z__1.i = 
		    z__2.r * du[i__3].i + z__2.i * du[i__3].r;
#line 206 "cgttrf.f"
	    du[i__2].r = z__1.r, du[i__2].i = z__1.i;
#line 207 "cgttrf.f"
	    ipiv[i__] = i__ + 1;
#line 208 "cgttrf.f"
	}
#line 209 "cgttrf.f"
/* L30: */
#line 209 "cgttrf.f"
    }
#line 210 "cgttrf.f"
    if (*n > 1) {
#line 211 "cgttrf.f"
	i__ = *n - 1;
#line 212 "cgttrf.f"
	i__1 = i__;
#line 212 "cgttrf.f"
	i__2 = i__;
#line 212 "cgttrf.f"
	if ((d__1 = d__[i__1].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) >= (d__3 = dl[i__2].r, abs(d__3)) + (d__4 = d_imag(&dl[
		i__]), abs(d__4))) {
#line 213 "cgttrf.f"
	    i__1 = i__;
#line 213 "cgttrf.f"
	    if ((d__1 = d__[i__1].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), 
		    abs(d__2)) != 0.) {
#line 214 "cgttrf.f"
		z_div(&z__1, &dl[i__], &d__[i__]);
#line 214 "cgttrf.f"
		fact.r = z__1.r, fact.i = z__1.i;
#line 215 "cgttrf.f"
		i__1 = i__;
#line 215 "cgttrf.f"
		dl[i__1].r = fact.r, dl[i__1].i = fact.i;
#line 216 "cgttrf.f"
		i__1 = i__ + 1;
#line 216 "cgttrf.f"
		i__2 = i__ + 1;
#line 216 "cgttrf.f"
		i__3 = i__;
#line 216 "cgttrf.f"
		z__2.r = fact.r * du[i__3].r - fact.i * du[i__3].i, z__2.i = 
			fact.r * du[i__3].i + fact.i * du[i__3].r;
#line 216 "cgttrf.f"
		z__1.r = d__[i__2].r - z__2.r, z__1.i = d__[i__2].i - z__2.i;
#line 216 "cgttrf.f"
		d__[i__1].r = z__1.r, d__[i__1].i = z__1.i;
#line 217 "cgttrf.f"
	    }
#line 218 "cgttrf.f"
	} else {
#line 219 "cgttrf.f"
	    z_div(&z__1, &d__[i__], &dl[i__]);
#line 219 "cgttrf.f"
	    fact.r = z__1.r, fact.i = z__1.i;
#line 220 "cgttrf.f"
	    i__1 = i__;
#line 220 "cgttrf.f"
	    i__2 = i__;
#line 220 "cgttrf.f"
	    d__[i__1].r = dl[i__2].r, d__[i__1].i = dl[i__2].i;
#line 221 "cgttrf.f"
	    i__1 = i__;
#line 221 "cgttrf.f"
	    dl[i__1].r = fact.r, dl[i__1].i = fact.i;
#line 222 "cgttrf.f"
	    i__1 = i__;
#line 222 "cgttrf.f"
	    temp.r = du[i__1].r, temp.i = du[i__1].i;
#line 223 "cgttrf.f"
	    i__1 = i__;
#line 223 "cgttrf.f"
	    i__2 = i__ + 1;
#line 223 "cgttrf.f"
	    du[i__1].r = d__[i__2].r, du[i__1].i = d__[i__2].i;
#line 224 "cgttrf.f"
	    i__1 = i__ + 1;
#line 224 "cgttrf.f"
	    i__2 = i__ + 1;
#line 224 "cgttrf.f"
	    z__2.r = fact.r * d__[i__2].r - fact.i * d__[i__2].i, z__2.i = 
		    fact.r * d__[i__2].i + fact.i * d__[i__2].r;
#line 224 "cgttrf.f"
	    z__1.r = temp.r - z__2.r, z__1.i = temp.i - z__2.i;
#line 224 "cgttrf.f"
	    d__[i__1].r = z__1.r, d__[i__1].i = z__1.i;
#line 225 "cgttrf.f"
	    ipiv[i__] = i__ + 1;
#line 226 "cgttrf.f"
	}
#line 227 "cgttrf.f"
    }

/*     Check for a zero on the diagonal of U. */

#line 231 "cgttrf.f"
    i__1 = *n;
#line 231 "cgttrf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 232 "cgttrf.f"
	i__2 = i__;
#line 232 "cgttrf.f"
	if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) == 0.) {
#line 233 "cgttrf.f"
	    *info = i__;
#line 234 "cgttrf.f"
	    goto L50;
#line 235 "cgttrf.f"
	}
#line 236 "cgttrf.f"
/* L40: */
#line 236 "cgttrf.f"
    }
#line 237 "cgttrf.f"
L50:

#line 239 "cgttrf.f"
    return 0;

/*     End of CGTTRF */

} /* cgttrf_ */

