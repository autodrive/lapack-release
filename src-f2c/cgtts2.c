#line 1 "cgtts2.f"
/* cgtts2.f -- translated by f2c (version 20100827).
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

#line 1 "cgtts2.f"
/* > \brief \b CGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization
 computed by sgttrf. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGTTS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtts2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtts2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtts2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ITRANS, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGTTS2 solves one of the systems of equations */
/* >    A * X = B,  A**T * X = B,  or  A**H * X = B, */
/* > with a tridiagonal matrix A using the LU factorization computed */
/* > by CGTTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITRANS */
/* > \verbatim */
/* >          ITRANS is INTEGER */
/* >          Specifies the form of the system of equations. */
/* >          = 0:  A * X = B     (No transpose) */
/* >          = 1:  A**T * X = B  (Transpose) */
/* >          = 2:  A**H * X = B  (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is COMPLEX array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is COMPLEX array, dimension (N-1) */
/* >          The (n-1) elements of the first super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is COMPLEX array, dimension (N-2) */
/* >          The (n-2) elements of the second super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= n, row i of the matrix was */
/* >          interchanged with row IPIV(i).  IPIV(i) will always be either */
/* >          i or i+1; IPIV(i) = i indicates a row interchange was not */
/* >          required. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
/* >          On entry, the matrix of right hand side vectors B. */
/* >          On exit, B is overwritten by the solution vectors X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgtts2_(integer *itrans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublecomplex temp;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 157 "cgtts2.f"
    /* Parameter adjustments */
#line 157 "cgtts2.f"
    --dl;
#line 157 "cgtts2.f"
    --d__;
#line 157 "cgtts2.f"
    --du;
#line 157 "cgtts2.f"
    --du2;
#line 157 "cgtts2.f"
    --ipiv;
#line 157 "cgtts2.f"
    b_dim1 = *ldb;
#line 157 "cgtts2.f"
    b_offset = 1 + b_dim1;
#line 157 "cgtts2.f"
    b -= b_offset;
#line 157 "cgtts2.f"

#line 157 "cgtts2.f"
    /* Function Body */
#line 157 "cgtts2.f"
    if (*n == 0 || *nrhs == 0) {
#line 157 "cgtts2.f"
	return 0;
#line 157 "cgtts2.f"
    }

#line 160 "cgtts2.f"
    if (*itrans == 0) {

/*        Solve A*X = B using the LU factorization of A, */
/*        overwriting each right hand side vector with its solution. */

#line 165 "cgtts2.f"
	if (*nrhs <= 1) {
#line 166 "cgtts2.f"
	    j = 1;
#line 167 "cgtts2.f"
L10:

/*           Solve L*x = b. */

#line 171 "cgtts2.f"
	    i__1 = *n - 1;
#line 171 "cgtts2.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 172 "cgtts2.f"
		if (ipiv[i__] == i__) {
#line 173 "cgtts2.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 173 "cgtts2.f"
		    i__3 = i__ + 1 + j * b_dim1;
#line 173 "cgtts2.f"
		    i__4 = i__;
#line 173 "cgtts2.f"
		    i__5 = i__ + j * b_dim1;
#line 173 "cgtts2.f"
		    z__2.r = dl[i__4].r * b[i__5].r - dl[i__4].i * b[i__5].i, 
			    z__2.i = dl[i__4].r * b[i__5].i + dl[i__4].i * b[
			    i__5].r;
#line 173 "cgtts2.f"
		    z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
#line 173 "cgtts2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 174 "cgtts2.f"
		} else {
#line 175 "cgtts2.f"
		    i__2 = i__ + j * b_dim1;
#line 175 "cgtts2.f"
		    temp.r = b[i__2].r, temp.i = b[i__2].i;
#line 176 "cgtts2.f"
		    i__2 = i__ + j * b_dim1;
#line 176 "cgtts2.f"
		    i__3 = i__ + 1 + j * b_dim1;
#line 176 "cgtts2.f"
		    b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
#line 177 "cgtts2.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 177 "cgtts2.f"
		    i__3 = i__;
#line 177 "cgtts2.f"
		    i__4 = i__ + j * b_dim1;
#line 177 "cgtts2.f"
		    z__2.r = dl[i__3].r * b[i__4].r - dl[i__3].i * b[i__4].i, 
			    z__2.i = dl[i__3].r * b[i__4].i + dl[i__3].i * b[
			    i__4].r;
#line 177 "cgtts2.f"
		    z__1.r = temp.r - z__2.r, z__1.i = temp.i - z__2.i;
#line 177 "cgtts2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 178 "cgtts2.f"
		}
#line 179 "cgtts2.f"
/* L20: */
#line 179 "cgtts2.f"
	    }

/*           Solve U*x = b. */

#line 183 "cgtts2.f"
	    i__1 = *n + j * b_dim1;
#line 183 "cgtts2.f"
	    z_div(&z__1, &b[*n + j * b_dim1], &d__[*n]);
#line 183 "cgtts2.f"
	    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 184 "cgtts2.f"
	    if (*n > 1) {
#line 184 "cgtts2.f"
		i__1 = *n - 1 + j * b_dim1;
#line 184 "cgtts2.f"
		i__2 = *n - 1 + j * b_dim1;
#line 184 "cgtts2.f"
		i__3 = *n - 1;
#line 184 "cgtts2.f"
		i__4 = *n + j * b_dim1;
#line 184 "cgtts2.f"
		z__3.r = du[i__3].r * b[i__4].r - du[i__3].i * b[i__4].i, 
			z__3.i = du[i__3].r * b[i__4].i + du[i__3].i * b[i__4]
			.r;
#line 184 "cgtts2.f"
		z__2.r = b[i__2].r - z__3.r, z__2.i = b[i__2].i - z__3.i;
#line 184 "cgtts2.f"
		z_div(&z__1, &z__2, &d__[*n - 1]);
#line 184 "cgtts2.f"
		b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 184 "cgtts2.f"
	    }
#line 187 "cgtts2.f"
	    for (i__ = *n - 2; i__ >= 1; --i__) {
#line 188 "cgtts2.f"
		i__1 = i__ + j * b_dim1;
#line 188 "cgtts2.f"
		i__2 = i__ + j * b_dim1;
#line 188 "cgtts2.f"
		i__3 = i__;
#line 188 "cgtts2.f"
		i__4 = i__ + 1 + j * b_dim1;
#line 188 "cgtts2.f"
		z__4.r = du[i__3].r * b[i__4].r - du[i__3].i * b[i__4].i, 
			z__4.i = du[i__3].r * b[i__4].i + du[i__3].i * b[i__4]
			.r;
#line 188 "cgtts2.f"
		z__3.r = b[i__2].r - z__4.r, z__3.i = b[i__2].i - z__4.i;
#line 188 "cgtts2.f"
		i__5 = i__;
#line 188 "cgtts2.f"
		i__6 = i__ + 2 + j * b_dim1;
#line 188 "cgtts2.f"
		z__5.r = du2[i__5].r * b[i__6].r - du2[i__5].i * b[i__6].i, 
			z__5.i = du2[i__5].r * b[i__6].i + du2[i__5].i * b[
			i__6].r;
#line 188 "cgtts2.f"
		z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
#line 188 "cgtts2.f"
		z_div(&z__1, &z__2, &d__[i__]);
#line 188 "cgtts2.f"
		b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 190 "cgtts2.f"
/* L30: */
#line 190 "cgtts2.f"
	    }
#line 191 "cgtts2.f"
	    if (j < *nrhs) {
#line 192 "cgtts2.f"
		++j;
#line 193 "cgtts2.f"
		goto L10;
#line 194 "cgtts2.f"
	    }
#line 195 "cgtts2.f"
	} else {
#line 196 "cgtts2.f"
	    i__1 = *nrhs;
#line 196 "cgtts2.f"
	    for (j = 1; j <= i__1; ++j) {

/*           Solve L*x = b. */

#line 200 "cgtts2.f"
		i__2 = *n - 1;
#line 200 "cgtts2.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 201 "cgtts2.f"
		    if (ipiv[i__] == i__) {
#line 202 "cgtts2.f"
			i__3 = i__ + 1 + j * b_dim1;
#line 202 "cgtts2.f"
			i__4 = i__ + 1 + j * b_dim1;
#line 202 "cgtts2.f"
			i__5 = i__;
#line 202 "cgtts2.f"
			i__6 = i__ + j * b_dim1;
#line 202 "cgtts2.f"
			z__2.r = dl[i__5].r * b[i__6].r - dl[i__5].i * b[i__6]
				.i, z__2.i = dl[i__5].r * b[i__6].i + dl[i__5]
				.i * b[i__6].r;
#line 202 "cgtts2.f"
			z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4].i - 
				z__2.i;
#line 202 "cgtts2.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 203 "cgtts2.f"
		    } else {
#line 204 "cgtts2.f"
			i__3 = i__ + j * b_dim1;
#line 204 "cgtts2.f"
			temp.r = b[i__3].r, temp.i = b[i__3].i;
#line 205 "cgtts2.f"
			i__3 = i__ + j * b_dim1;
#line 205 "cgtts2.f"
			i__4 = i__ + 1 + j * b_dim1;
#line 205 "cgtts2.f"
			b[i__3].r = b[i__4].r, b[i__3].i = b[i__4].i;
#line 206 "cgtts2.f"
			i__3 = i__ + 1 + j * b_dim1;
#line 206 "cgtts2.f"
			i__4 = i__;
#line 206 "cgtts2.f"
			i__5 = i__ + j * b_dim1;
#line 206 "cgtts2.f"
			z__2.r = dl[i__4].r * b[i__5].r - dl[i__4].i * b[i__5]
				.i, z__2.i = dl[i__4].r * b[i__5].i + dl[i__4]
				.i * b[i__5].r;
#line 206 "cgtts2.f"
			z__1.r = temp.r - z__2.r, z__1.i = temp.i - z__2.i;
#line 206 "cgtts2.f"
			b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 207 "cgtts2.f"
		    }
#line 208 "cgtts2.f"
/* L40: */
#line 208 "cgtts2.f"
		}

/*           Solve U*x = b. */

#line 212 "cgtts2.f"
		i__2 = *n + j * b_dim1;
#line 212 "cgtts2.f"
		z_div(&z__1, &b[*n + j * b_dim1], &d__[*n]);
#line 212 "cgtts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 213 "cgtts2.f"
		if (*n > 1) {
#line 213 "cgtts2.f"
		    i__2 = *n - 1 + j * b_dim1;
#line 213 "cgtts2.f"
		    i__3 = *n - 1 + j * b_dim1;
#line 213 "cgtts2.f"
		    i__4 = *n - 1;
#line 213 "cgtts2.f"
		    i__5 = *n + j * b_dim1;
#line 213 "cgtts2.f"
		    z__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, 
			    z__3.i = du[i__4].r * b[i__5].i + du[i__4].i * b[
			    i__5].r;
#line 213 "cgtts2.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 213 "cgtts2.f"
		    z_div(&z__1, &z__2, &d__[*n - 1]);
#line 213 "cgtts2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 213 "cgtts2.f"
		}
#line 216 "cgtts2.f"
		for (i__ = *n - 2; i__ >= 1; --i__) {
#line 217 "cgtts2.f"
		    i__2 = i__ + j * b_dim1;
#line 217 "cgtts2.f"
		    i__3 = i__ + j * b_dim1;
#line 217 "cgtts2.f"
		    i__4 = i__;
#line 217 "cgtts2.f"
		    i__5 = i__ + 1 + j * b_dim1;
#line 217 "cgtts2.f"
		    z__4.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, 
			    z__4.i = du[i__4].r * b[i__5].i + du[i__4].i * b[
			    i__5].r;
#line 217 "cgtts2.f"
		    z__3.r = b[i__3].r - z__4.r, z__3.i = b[i__3].i - z__4.i;
#line 217 "cgtts2.f"
		    i__6 = i__;
#line 217 "cgtts2.f"
		    i__7 = i__ + 2 + j * b_dim1;
#line 217 "cgtts2.f"
		    z__5.r = du2[i__6].r * b[i__7].r - du2[i__6].i * b[i__7]
			    .i, z__5.i = du2[i__6].r * b[i__7].i + du2[i__6]
			    .i * b[i__7].r;
#line 217 "cgtts2.f"
		    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
#line 217 "cgtts2.f"
		    z_div(&z__1, &z__2, &d__[i__]);
#line 217 "cgtts2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 219 "cgtts2.f"
/* L50: */
#line 219 "cgtts2.f"
		}
#line 220 "cgtts2.f"
/* L60: */
#line 220 "cgtts2.f"
	    }
#line 221 "cgtts2.f"
	}
#line 222 "cgtts2.f"
    } else if (*itrans == 1) {

/*        Solve A**T * X = B. */

#line 226 "cgtts2.f"
	if (*nrhs <= 1) {
#line 227 "cgtts2.f"
	    j = 1;
#line 228 "cgtts2.f"
L70:

/*           Solve U**T * x = b. */

#line 232 "cgtts2.f"
	    i__1 = j * b_dim1 + 1;
#line 232 "cgtts2.f"
	    z_div(&z__1, &b[j * b_dim1 + 1], &d__[1]);
#line 232 "cgtts2.f"
	    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 233 "cgtts2.f"
	    if (*n > 1) {
#line 233 "cgtts2.f"
		i__1 = j * b_dim1 + 2;
#line 233 "cgtts2.f"
		i__2 = j * b_dim1 + 2;
#line 233 "cgtts2.f"
		i__3 = j * b_dim1 + 1;
#line 233 "cgtts2.f"
		z__3.r = du[1].r * b[i__3].r - du[1].i * b[i__3].i, z__3.i = 
			du[1].r * b[i__3].i + du[1].i * b[i__3].r;
#line 233 "cgtts2.f"
		z__2.r = b[i__2].r - z__3.r, z__2.i = b[i__2].i - z__3.i;
#line 233 "cgtts2.f"
		z_div(&z__1, &z__2, &d__[2]);
#line 233 "cgtts2.f"
		b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 233 "cgtts2.f"
	    }
#line 235 "cgtts2.f"
	    i__1 = *n;
#line 235 "cgtts2.f"
	    for (i__ = 3; i__ <= i__1; ++i__) {
#line 236 "cgtts2.f"
		i__2 = i__ + j * b_dim1;
#line 236 "cgtts2.f"
		i__3 = i__ + j * b_dim1;
#line 236 "cgtts2.f"
		i__4 = i__ - 1;
#line 236 "cgtts2.f"
		i__5 = i__ - 1 + j * b_dim1;
#line 236 "cgtts2.f"
		z__4.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, 
			z__4.i = du[i__4].r * b[i__5].i + du[i__4].i * b[i__5]
			.r;
#line 236 "cgtts2.f"
		z__3.r = b[i__3].r - z__4.r, z__3.i = b[i__3].i - z__4.i;
#line 236 "cgtts2.f"
		i__6 = i__ - 2;
#line 236 "cgtts2.f"
		i__7 = i__ - 2 + j * b_dim1;
#line 236 "cgtts2.f"
		z__5.r = du2[i__6].r * b[i__7].r - du2[i__6].i * b[i__7].i, 
			z__5.i = du2[i__6].r * b[i__7].i + du2[i__6].i * b[
			i__7].r;
#line 236 "cgtts2.f"
		z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
#line 236 "cgtts2.f"
		z_div(&z__1, &z__2, &d__[i__]);
#line 236 "cgtts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 238 "cgtts2.f"
/* L80: */
#line 238 "cgtts2.f"
	    }

/*           Solve L**T * x = b. */

#line 242 "cgtts2.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 243 "cgtts2.f"
		if (ipiv[i__] == i__) {
#line 244 "cgtts2.f"
		    i__1 = i__ + j * b_dim1;
#line 244 "cgtts2.f"
		    i__2 = i__ + j * b_dim1;
#line 244 "cgtts2.f"
		    i__3 = i__;
#line 244 "cgtts2.f"
		    i__4 = i__ + 1 + j * b_dim1;
#line 244 "cgtts2.f"
		    z__2.r = dl[i__3].r * b[i__4].r - dl[i__3].i * b[i__4].i, 
			    z__2.i = dl[i__3].r * b[i__4].i + dl[i__3].i * b[
			    i__4].r;
#line 244 "cgtts2.f"
		    z__1.r = b[i__2].r - z__2.r, z__1.i = b[i__2].i - z__2.i;
#line 244 "cgtts2.f"
		    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 245 "cgtts2.f"
		} else {
#line 246 "cgtts2.f"
		    i__1 = i__ + 1 + j * b_dim1;
#line 246 "cgtts2.f"
		    temp.r = b[i__1].r, temp.i = b[i__1].i;
#line 247 "cgtts2.f"
		    i__1 = i__ + 1 + j * b_dim1;
#line 247 "cgtts2.f"
		    i__2 = i__ + j * b_dim1;
#line 247 "cgtts2.f"
		    i__3 = i__;
#line 247 "cgtts2.f"
		    z__2.r = dl[i__3].r * temp.r - dl[i__3].i * temp.i, 
			    z__2.i = dl[i__3].r * temp.i + dl[i__3].i * 
			    temp.r;
#line 247 "cgtts2.f"
		    z__1.r = b[i__2].r - z__2.r, z__1.i = b[i__2].i - z__2.i;
#line 247 "cgtts2.f"
		    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 248 "cgtts2.f"
		    i__1 = i__ + j * b_dim1;
#line 248 "cgtts2.f"
		    b[i__1].r = temp.r, b[i__1].i = temp.i;
#line 249 "cgtts2.f"
		}
#line 250 "cgtts2.f"
/* L90: */
#line 250 "cgtts2.f"
	    }
#line 251 "cgtts2.f"
	    if (j < *nrhs) {
#line 252 "cgtts2.f"
		++j;
#line 253 "cgtts2.f"
		goto L70;
#line 254 "cgtts2.f"
	    }
#line 255 "cgtts2.f"
	} else {
#line 256 "cgtts2.f"
	    i__1 = *nrhs;
#line 256 "cgtts2.f"
	    for (j = 1; j <= i__1; ++j) {

/*           Solve U**T * x = b. */

#line 260 "cgtts2.f"
		i__2 = j * b_dim1 + 1;
#line 260 "cgtts2.f"
		z_div(&z__1, &b[j * b_dim1 + 1], &d__[1]);
#line 260 "cgtts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 261 "cgtts2.f"
		if (*n > 1) {
#line 261 "cgtts2.f"
		    i__2 = j * b_dim1 + 2;
#line 261 "cgtts2.f"
		    i__3 = j * b_dim1 + 2;
#line 261 "cgtts2.f"
		    i__4 = j * b_dim1 + 1;
#line 261 "cgtts2.f"
		    z__3.r = du[1].r * b[i__4].r - du[1].i * b[i__4].i, 
			    z__3.i = du[1].r * b[i__4].i + du[1].i * b[i__4]
			    .r;
#line 261 "cgtts2.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 261 "cgtts2.f"
		    z_div(&z__1, &z__2, &d__[2]);
#line 261 "cgtts2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 261 "cgtts2.f"
		}
#line 263 "cgtts2.f"
		i__2 = *n;
#line 263 "cgtts2.f"
		for (i__ = 3; i__ <= i__2; ++i__) {
#line 264 "cgtts2.f"
		    i__3 = i__ + j * b_dim1;
#line 264 "cgtts2.f"
		    i__4 = i__ + j * b_dim1;
#line 264 "cgtts2.f"
		    i__5 = i__ - 1;
#line 264 "cgtts2.f"
		    i__6 = i__ - 1 + j * b_dim1;
#line 264 "cgtts2.f"
		    z__4.r = du[i__5].r * b[i__6].r - du[i__5].i * b[i__6].i, 
			    z__4.i = du[i__5].r * b[i__6].i + du[i__5].i * b[
			    i__6].r;
#line 264 "cgtts2.f"
		    z__3.r = b[i__4].r - z__4.r, z__3.i = b[i__4].i - z__4.i;
#line 264 "cgtts2.f"
		    i__7 = i__ - 2;
#line 264 "cgtts2.f"
		    i__8 = i__ - 2 + j * b_dim1;
#line 264 "cgtts2.f"
		    z__5.r = du2[i__7].r * b[i__8].r - du2[i__7].i * b[i__8]
			    .i, z__5.i = du2[i__7].r * b[i__8].i + du2[i__7]
			    .i * b[i__8].r;
#line 264 "cgtts2.f"
		    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
#line 264 "cgtts2.f"
		    z_div(&z__1, &z__2, &d__[i__]);
#line 264 "cgtts2.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 266 "cgtts2.f"
/* L100: */
#line 266 "cgtts2.f"
		}

/*           Solve L**T * x = b. */

#line 270 "cgtts2.f"
		for (i__ = *n - 1; i__ >= 1; --i__) {
#line 271 "cgtts2.f"
		    if (ipiv[i__] == i__) {
#line 272 "cgtts2.f"
			i__2 = i__ + j * b_dim1;
#line 272 "cgtts2.f"
			i__3 = i__ + j * b_dim1;
#line 272 "cgtts2.f"
			i__4 = i__;
#line 272 "cgtts2.f"
			i__5 = i__ + 1 + j * b_dim1;
#line 272 "cgtts2.f"
			z__2.r = dl[i__4].r * b[i__5].r - dl[i__4].i * b[i__5]
				.i, z__2.i = dl[i__4].r * b[i__5].i + dl[i__4]
				.i * b[i__5].r;
#line 272 "cgtts2.f"
			z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - 
				z__2.i;
#line 272 "cgtts2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 273 "cgtts2.f"
		    } else {
#line 274 "cgtts2.f"
			i__2 = i__ + 1 + j * b_dim1;
#line 274 "cgtts2.f"
			temp.r = b[i__2].r, temp.i = b[i__2].i;
#line 275 "cgtts2.f"
			i__2 = i__ + 1 + j * b_dim1;
#line 275 "cgtts2.f"
			i__3 = i__ + j * b_dim1;
#line 275 "cgtts2.f"
			i__4 = i__;
#line 275 "cgtts2.f"
			z__2.r = dl[i__4].r * temp.r - dl[i__4].i * temp.i, 
				z__2.i = dl[i__4].r * temp.i + dl[i__4].i * 
				temp.r;
#line 275 "cgtts2.f"
			z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - 
				z__2.i;
#line 275 "cgtts2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 276 "cgtts2.f"
			i__2 = i__ + j * b_dim1;
#line 276 "cgtts2.f"
			b[i__2].r = temp.r, b[i__2].i = temp.i;
#line 277 "cgtts2.f"
		    }
#line 278 "cgtts2.f"
/* L110: */
#line 278 "cgtts2.f"
		}
#line 279 "cgtts2.f"
/* L120: */
#line 279 "cgtts2.f"
	    }
#line 280 "cgtts2.f"
	}
#line 281 "cgtts2.f"
    } else {

/*        Solve A**H * X = B. */

#line 285 "cgtts2.f"
	if (*nrhs <= 1) {
#line 286 "cgtts2.f"
	    j = 1;
#line 287 "cgtts2.f"
L130:

/*           Solve U**H * x = b. */

#line 291 "cgtts2.f"
	    i__1 = j * b_dim1 + 1;
#line 291 "cgtts2.f"
	    d_cnjg(&z__2, &d__[1]);
#line 291 "cgtts2.f"
	    z_div(&z__1, &b[j * b_dim1 + 1], &z__2);
#line 291 "cgtts2.f"
	    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 292 "cgtts2.f"
	    if (*n > 1) {
#line 292 "cgtts2.f"
		i__1 = j * b_dim1 + 2;
#line 292 "cgtts2.f"
		i__2 = j * b_dim1 + 2;
#line 292 "cgtts2.f"
		d_cnjg(&z__4, &du[1]);
#line 292 "cgtts2.f"
		i__3 = j * b_dim1 + 1;
#line 292 "cgtts2.f"
		z__3.r = z__4.r * b[i__3].r - z__4.i * b[i__3].i, z__3.i = 
			z__4.r * b[i__3].i + z__4.i * b[i__3].r;
#line 292 "cgtts2.f"
		z__2.r = b[i__2].r - z__3.r, z__2.i = b[i__2].i - z__3.i;
#line 292 "cgtts2.f"
		d_cnjg(&z__5, &d__[2]);
#line 292 "cgtts2.f"
		z_div(&z__1, &z__2, &z__5);
#line 292 "cgtts2.f"
		b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 292 "cgtts2.f"
	    }
#line 295 "cgtts2.f"
	    i__1 = *n;
#line 295 "cgtts2.f"
	    for (i__ = 3; i__ <= i__1; ++i__) {
#line 296 "cgtts2.f"
		i__2 = i__ + j * b_dim1;
#line 296 "cgtts2.f"
		i__3 = i__ + j * b_dim1;
#line 296 "cgtts2.f"
		d_cnjg(&z__5, &du[i__ - 1]);
#line 296 "cgtts2.f"
		i__4 = i__ - 1 + j * b_dim1;
#line 296 "cgtts2.f"
		z__4.r = z__5.r * b[i__4].r - z__5.i * b[i__4].i, z__4.i = 
			z__5.r * b[i__4].i + z__5.i * b[i__4].r;
#line 296 "cgtts2.f"
		z__3.r = b[i__3].r - z__4.r, z__3.i = b[i__3].i - z__4.i;
#line 296 "cgtts2.f"
		d_cnjg(&z__7, &du2[i__ - 2]);
#line 296 "cgtts2.f"
		i__5 = i__ - 2 + j * b_dim1;
#line 296 "cgtts2.f"
		z__6.r = z__7.r * b[i__5].r - z__7.i * b[i__5].i, z__6.i = 
			z__7.r * b[i__5].i + z__7.i * b[i__5].r;
#line 296 "cgtts2.f"
		z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 296 "cgtts2.f"
		d_cnjg(&z__8, &d__[i__]);
#line 296 "cgtts2.f"
		z_div(&z__1, &z__2, &z__8);
#line 296 "cgtts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 299 "cgtts2.f"
/* L140: */
#line 299 "cgtts2.f"
	    }

/*           Solve L**H * x = b. */

#line 303 "cgtts2.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 304 "cgtts2.f"
		if (ipiv[i__] == i__) {
#line 305 "cgtts2.f"
		    i__1 = i__ + j * b_dim1;
#line 305 "cgtts2.f"
		    i__2 = i__ + j * b_dim1;
#line 305 "cgtts2.f"
		    d_cnjg(&z__3, &dl[i__]);
#line 305 "cgtts2.f"
		    i__3 = i__ + 1 + j * b_dim1;
#line 305 "cgtts2.f"
		    z__2.r = z__3.r * b[i__3].r - z__3.i * b[i__3].i, z__2.i =
			     z__3.r * b[i__3].i + z__3.i * b[i__3].r;
#line 305 "cgtts2.f"
		    z__1.r = b[i__2].r - z__2.r, z__1.i = b[i__2].i - z__2.i;
#line 305 "cgtts2.f"
		    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 306 "cgtts2.f"
		} else {
#line 307 "cgtts2.f"
		    i__1 = i__ + 1 + j * b_dim1;
#line 307 "cgtts2.f"
		    temp.r = b[i__1].r, temp.i = b[i__1].i;
#line 308 "cgtts2.f"
		    i__1 = i__ + 1 + j * b_dim1;
#line 308 "cgtts2.f"
		    i__2 = i__ + j * b_dim1;
#line 308 "cgtts2.f"
		    d_cnjg(&z__3, &dl[i__]);
#line 308 "cgtts2.f"
		    z__2.r = z__3.r * temp.r - z__3.i * temp.i, z__2.i = 
			    z__3.r * temp.i + z__3.i * temp.r;
#line 308 "cgtts2.f"
		    z__1.r = b[i__2].r - z__2.r, z__1.i = b[i__2].i - z__2.i;
#line 308 "cgtts2.f"
		    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 309 "cgtts2.f"
		    i__1 = i__ + j * b_dim1;
#line 309 "cgtts2.f"
		    b[i__1].r = temp.r, b[i__1].i = temp.i;
#line 310 "cgtts2.f"
		}
#line 311 "cgtts2.f"
/* L150: */
#line 311 "cgtts2.f"
	    }
#line 312 "cgtts2.f"
	    if (j < *nrhs) {
#line 313 "cgtts2.f"
		++j;
#line 314 "cgtts2.f"
		goto L130;
#line 315 "cgtts2.f"
	    }
#line 316 "cgtts2.f"
	} else {
#line 317 "cgtts2.f"
	    i__1 = *nrhs;
#line 317 "cgtts2.f"
	    for (j = 1; j <= i__1; ++j) {

/*           Solve U**H * x = b. */

#line 321 "cgtts2.f"
		i__2 = j * b_dim1 + 1;
#line 321 "cgtts2.f"
		d_cnjg(&z__2, &d__[1]);
#line 321 "cgtts2.f"
		z_div(&z__1, &b[j * b_dim1 + 1], &z__2);
#line 321 "cgtts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 322 "cgtts2.f"
		if (*n > 1) {
#line 322 "cgtts2.f"
		    i__2 = j * b_dim1 + 2;
#line 322 "cgtts2.f"
		    i__3 = j * b_dim1 + 2;
#line 322 "cgtts2.f"
		    d_cnjg(&z__4, &du[1]);
#line 322 "cgtts2.f"
		    i__4 = j * b_dim1 + 1;
#line 322 "cgtts2.f"
		    z__3.r = z__4.r * b[i__4].r - z__4.i * b[i__4].i, z__3.i =
			     z__4.r * b[i__4].i + z__4.i * b[i__4].r;
#line 322 "cgtts2.f"
		    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 322 "cgtts2.f"
		    d_cnjg(&z__5, &d__[2]);
#line 322 "cgtts2.f"
		    z_div(&z__1, &z__2, &z__5);
#line 322 "cgtts2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 322 "cgtts2.f"
		}
#line 325 "cgtts2.f"
		i__2 = *n;
#line 325 "cgtts2.f"
		for (i__ = 3; i__ <= i__2; ++i__) {
#line 326 "cgtts2.f"
		    i__3 = i__ + j * b_dim1;
#line 326 "cgtts2.f"
		    i__4 = i__ + j * b_dim1;
#line 326 "cgtts2.f"
		    d_cnjg(&z__5, &du[i__ - 1]);
#line 326 "cgtts2.f"
		    i__5 = i__ - 1 + j * b_dim1;
#line 326 "cgtts2.f"
		    z__4.r = z__5.r * b[i__5].r - z__5.i * b[i__5].i, z__4.i =
			     z__5.r * b[i__5].i + z__5.i * b[i__5].r;
#line 326 "cgtts2.f"
		    z__3.r = b[i__4].r - z__4.r, z__3.i = b[i__4].i - z__4.i;
#line 326 "cgtts2.f"
		    d_cnjg(&z__7, &du2[i__ - 2]);
#line 326 "cgtts2.f"
		    i__6 = i__ - 2 + j * b_dim1;
#line 326 "cgtts2.f"
		    z__6.r = z__7.r * b[i__6].r - z__7.i * b[i__6].i, z__6.i =
			     z__7.r * b[i__6].i + z__7.i * b[i__6].r;
#line 326 "cgtts2.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 326 "cgtts2.f"
		    d_cnjg(&z__8, &d__[i__]);
#line 326 "cgtts2.f"
		    z_div(&z__1, &z__2, &z__8);
#line 326 "cgtts2.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 329 "cgtts2.f"
/* L160: */
#line 329 "cgtts2.f"
		}

/*           Solve L**H * x = b. */

#line 333 "cgtts2.f"
		for (i__ = *n - 1; i__ >= 1; --i__) {
#line 334 "cgtts2.f"
		    if (ipiv[i__] == i__) {
#line 335 "cgtts2.f"
			i__2 = i__ + j * b_dim1;
#line 335 "cgtts2.f"
			i__3 = i__ + j * b_dim1;
#line 335 "cgtts2.f"
			d_cnjg(&z__3, &dl[i__]);
#line 335 "cgtts2.f"
			i__4 = i__ + 1 + j * b_dim1;
#line 335 "cgtts2.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 335 "cgtts2.f"
			z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - 
				z__2.i;
#line 335 "cgtts2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 337 "cgtts2.f"
		    } else {
#line 338 "cgtts2.f"
			i__2 = i__ + 1 + j * b_dim1;
#line 338 "cgtts2.f"
			temp.r = b[i__2].r, temp.i = b[i__2].i;
#line 339 "cgtts2.f"
			i__2 = i__ + 1 + j * b_dim1;
#line 339 "cgtts2.f"
			i__3 = i__ + j * b_dim1;
#line 339 "cgtts2.f"
			d_cnjg(&z__3, &dl[i__]);
#line 339 "cgtts2.f"
			z__2.r = z__3.r * temp.r - z__3.i * temp.i, z__2.i = 
				z__3.r * temp.i + z__3.i * temp.r;
#line 339 "cgtts2.f"
			z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - 
				z__2.i;
#line 339 "cgtts2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 340 "cgtts2.f"
			i__2 = i__ + j * b_dim1;
#line 340 "cgtts2.f"
			b[i__2].r = temp.r, b[i__2].i = temp.i;
#line 341 "cgtts2.f"
		    }
#line 342 "cgtts2.f"
/* L170: */
#line 342 "cgtts2.f"
		}
#line 343 "cgtts2.f"
/* L180: */
#line 343 "cgtts2.f"
	    }
#line 344 "cgtts2.f"
	}
#line 345 "cgtts2.f"
    }

/*     End of CGTTS2 */

#line 349 "cgtts2.f"
    return 0;
} /* cgtts2_ */

