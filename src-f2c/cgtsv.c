#line 1 "cgtsv.f"
/* cgtsv.f -- translated by f2c (version 20100827).
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

#line 1 "cgtsv.f"
/* > \brief <b> CGTSV computes the solution to system of linear equations A * X = B for GT matrices </b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGTSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtsv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtsv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtsv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGTSV( N, NRHS, DL, D, DU, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGTSV  solves the equation */
/* > */
/* >    A*X = B, */
/* > */
/* > where A is an N-by-N tridiagonal matrix, by Gaussian elimination with */
/* > partial pivoting. */
/* > */
/* > Note that the equation  A**T *X = B  may be solved by interchanging the */
/* > order of the arguments DU and DL. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* >          DL is COMPLEX array, dimension (N-1) */
/* >          On entry, DL must contain the (n-1) subdiagonal elements of */
/* >          A. */
/* >          On exit, DL is overwritten by the (n-2) elements of the */
/* >          second superdiagonal of the upper triangular matrix U from */
/* >          the LU factorization of A, in DL(1), ..., DL(n-2). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          On entry, D must contain the diagonal elements of A. */
/* >          On exit, D is overwritten by the n diagonal elements of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* >          DU is COMPLEX array, dimension (N-1) */
/* >          On entry, DU must contain the (n-1) superdiagonal elements */
/* >          of A. */
/* >          On exit, DU is overwritten by the (n-1) elements of the first */
/* >          superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
/* >          On entry, the N-by-NRHS right hand side matrix B. */
/* >          On exit, if INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, U(i,i) is exactly zero, and the solution */
/* >                has not been computed.  The factorization has not been */
/* >                completed unless i = N. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexGTsolve */

/*  ===================================================================== */
/* Subroutine */ int cgtsv_(integer *n, integer *nrhs, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *b, integer *ldb,
	 integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex temp, mult;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 163 "cgtsv.f"
    /* Parameter adjustments */
#line 163 "cgtsv.f"
    --dl;
#line 163 "cgtsv.f"
    --d__;
#line 163 "cgtsv.f"
    --du;
#line 163 "cgtsv.f"
    b_dim1 = *ldb;
#line 163 "cgtsv.f"
    b_offset = 1 + b_dim1;
#line 163 "cgtsv.f"
    b -= b_offset;
#line 163 "cgtsv.f"

#line 163 "cgtsv.f"
    /* Function Body */
#line 163 "cgtsv.f"
    *info = 0;
#line 164 "cgtsv.f"
    if (*n < 0) {
#line 165 "cgtsv.f"
	*info = -1;
#line 166 "cgtsv.f"
    } else if (*nrhs < 0) {
#line 167 "cgtsv.f"
	*info = -2;
#line 168 "cgtsv.f"
    } else if (*ldb < max(1,*n)) {
#line 169 "cgtsv.f"
	*info = -7;
#line 170 "cgtsv.f"
    }
#line 171 "cgtsv.f"
    if (*info != 0) {
#line 172 "cgtsv.f"
	i__1 = -(*info);
#line 172 "cgtsv.f"
	xerbla_("CGTSV ", &i__1, (ftnlen)6);
#line 173 "cgtsv.f"
	return 0;
#line 174 "cgtsv.f"
    }

#line 176 "cgtsv.f"
    if (*n == 0) {
#line 176 "cgtsv.f"
	return 0;
#line 176 "cgtsv.f"
    }

#line 179 "cgtsv.f"
    i__1 = *n - 1;
#line 179 "cgtsv.f"
    for (k = 1; k <= i__1; ++k) {
#line 180 "cgtsv.f"
	i__2 = k;
#line 180 "cgtsv.f"
	if (dl[i__2].r == 0. && dl[i__2].i == 0.) {

/*           Subdiagonal is zero, no elimination is required. */

#line 184 "cgtsv.f"
	    i__2 = k;
#line 184 "cgtsv.f"
	    if (d__[i__2].r == 0. && d__[i__2].i == 0.) {

/*              Diagonal is zero: set INFO = K and return; a unique */
/*              solution can not be found. */

#line 189 "cgtsv.f"
		*info = k;
#line 190 "cgtsv.f"
		return 0;
#line 191 "cgtsv.f"
	    }
#line 192 "cgtsv.f"
	} else /* if(complicated condition) */ {
#line 192 "cgtsv.f"
	    i__2 = k;
#line 192 "cgtsv.f"
	    i__3 = k;
#line 192 "cgtsv.f"
	    if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[k]), 
		    abs(d__2)) >= (d__3 = dl[i__3].r, abs(d__3)) + (d__4 = 
		    d_imag(&dl[k]), abs(d__4))) {

/*           No row interchange required */

#line 196 "cgtsv.f"
		z_div(&z__1, &dl[k], &d__[k]);
#line 196 "cgtsv.f"
		mult.r = z__1.r, mult.i = z__1.i;
#line 197 "cgtsv.f"
		i__2 = k + 1;
#line 197 "cgtsv.f"
		i__3 = k + 1;
#line 197 "cgtsv.f"
		i__4 = k;
#line 197 "cgtsv.f"
		z__2.r = mult.r * du[i__4].r - mult.i * du[i__4].i, z__2.i = 
			mult.r * du[i__4].i + mult.i * du[i__4].r;
#line 197 "cgtsv.f"
		z__1.r = d__[i__3].r - z__2.r, z__1.i = d__[i__3].i - z__2.i;
#line 197 "cgtsv.f"
		d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
#line 198 "cgtsv.f"
		i__2 = *nrhs;
#line 198 "cgtsv.f"
		for (j = 1; j <= i__2; ++j) {
#line 199 "cgtsv.f"
		    i__3 = k + 1 + j * b_dim1;
#line 199 "cgtsv.f"
		    i__4 = k + 1 + j * b_dim1;
#line 199 "cgtsv.f"
		    i__5 = k + j * b_dim1;
#line 199 "cgtsv.f"
		    z__2.r = mult.r * b[i__5].r - mult.i * b[i__5].i, z__2.i =
			     mult.r * b[i__5].i + mult.i * b[i__5].r;
#line 199 "cgtsv.f"
		    z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4].i - z__2.i;
#line 199 "cgtsv.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 200 "cgtsv.f"
/* L10: */
#line 200 "cgtsv.f"
		}
#line 201 "cgtsv.f"
		if (k < *n - 1) {
#line 201 "cgtsv.f"
		    i__2 = k;
#line 201 "cgtsv.f"
		    dl[i__2].r = 0., dl[i__2].i = 0.;
#line 201 "cgtsv.f"
		}
#line 203 "cgtsv.f"
	    } else {

/*           Interchange rows K and K+1 */

#line 207 "cgtsv.f"
		z_div(&z__1, &d__[k], &dl[k]);
#line 207 "cgtsv.f"
		mult.r = z__1.r, mult.i = z__1.i;
#line 208 "cgtsv.f"
		i__2 = k;
#line 208 "cgtsv.f"
		i__3 = k;
#line 208 "cgtsv.f"
		d__[i__2].r = dl[i__3].r, d__[i__2].i = dl[i__3].i;
#line 209 "cgtsv.f"
		i__2 = k + 1;
#line 209 "cgtsv.f"
		temp.r = d__[i__2].r, temp.i = d__[i__2].i;
#line 210 "cgtsv.f"
		i__2 = k + 1;
#line 210 "cgtsv.f"
		i__3 = k;
#line 210 "cgtsv.f"
		z__2.r = mult.r * temp.r - mult.i * temp.i, z__2.i = mult.r * 
			temp.i + mult.i * temp.r;
#line 210 "cgtsv.f"
		z__1.r = du[i__3].r - z__2.r, z__1.i = du[i__3].i - z__2.i;
#line 210 "cgtsv.f"
		d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
#line 211 "cgtsv.f"
		if (k < *n - 1) {
#line 212 "cgtsv.f"
		    i__2 = k;
#line 212 "cgtsv.f"
		    i__3 = k + 1;
#line 212 "cgtsv.f"
		    dl[i__2].r = du[i__3].r, dl[i__2].i = du[i__3].i;
#line 213 "cgtsv.f"
		    i__2 = k + 1;
#line 213 "cgtsv.f"
		    z__2.r = -mult.r, z__2.i = -mult.i;
#line 213 "cgtsv.f"
		    i__3 = k;
#line 213 "cgtsv.f"
		    z__1.r = z__2.r * dl[i__3].r - z__2.i * dl[i__3].i, 
			    z__1.i = z__2.r * dl[i__3].i + z__2.i * dl[i__3]
			    .r;
#line 213 "cgtsv.f"
		    du[i__2].r = z__1.r, du[i__2].i = z__1.i;
#line 214 "cgtsv.f"
		}
#line 215 "cgtsv.f"
		i__2 = k;
#line 215 "cgtsv.f"
		du[i__2].r = temp.r, du[i__2].i = temp.i;
#line 216 "cgtsv.f"
		i__2 = *nrhs;
#line 216 "cgtsv.f"
		for (j = 1; j <= i__2; ++j) {
#line 217 "cgtsv.f"
		    i__3 = k + j * b_dim1;
#line 217 "cgtsv.f"
		    temp.r = b[i__3].r, temp.i = b[i__3].i;
#line 218 "cgtsv.f"
		    i__3 = k + j * b_dim1;
#line 218 "cgtsv.f"
		    i__4 = k + 1 + j * b_dim1;
#line 218 "cgtsv.f"
		    b[i__3].r = b[i__4].r, b[i__3].i = b[i__4].i;
#line 219 "cgtsv.f"
		    i__3 = k + 1 + j * b_dim1;
#line 219 "cgtsv.f"
		    i__4 = k + 1 + j * b_dim1;
#line 219 "cgtsv.f"
		    z__2.r = mult.r * b[i__4].r - mult.i * b[i__4].i, z__2.i =
			     mult.r * b[i__4].i + mult.i * b[i__4].r;
#line 219 "cgtsv.f"
		    z__1.r = temp.r - z__2.r, z__1.i = temp.i - z__2.i;
#line 219 "cgtsv.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 220 "cgtsv.f"
/* L20: */
#line 220 "cgtsv.f"
		}
#line 221 "cgtsv.f"
	    }
#line 221 "cgtsv.f"
	}
#line 222 "cgtsv.f"
/* L30: */
#line 222 "cgtsv.f"
    }
#line 223 "cgtsv.f"
    i__1 = *n;
#line 223 "cgtsv.f"
    if (d__[i__1].r == 0. && d__[i__1].i == 0.) {
#line 224 "cgtsv.f"
	*info = *n;
#line 225 "cgtsv.f"
	return 0;
#line 226 "cgtsv.f"
    }

/*     Back solve with the matrix U from the factorization. */

#line 230 "cgtsv.f"
    i__1 = *nrhs;
#line 230 "cgtsv.f"
    for (j = 1; j <= i__1; ++j) {
#line 231 "cgtsv.f"
	i__2 = *n + j * b_dim1;
#line 231 "cgtsv.f"
	z_div(&z__1, &b[*n + j * b_dim1], &d__[*n]);
#line 231 "cgtsv.f"
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 232 "cgtsv.f"
	if (*n > 1) {
#line 232 "cgtsv.f"
	    i__2 = *n - 1 + j * b_dim1;
#line 232 "cgtsv.f"
	    i__3 = *n - 1 + j * b_dim1;
#line 232 "cgtsv.f"
	    i__4 = *n - 1;
#line 232 "cgtsv.f"
	    i__5 = *n + j * b_dim1;
#line 232 "cgtsv.f"
	    z__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, z__3.i =
		     du[i__4].r * b[i__5].i + du[i__4].i * b[i__5].r;
#line 232 "cgtsv.f"
	    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
#line 232 "cgtsv.f"
	    z_div(&z__1, &z__2, &d__[*n - 1]);
#line 232 "cgtsv.f"
	    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 232 "cgtsv.f"
	}
#line 234 "cgtsv.f"
	for (k = *n - 2; k >= 1; --k) {
#line 235 "cgtsv.f"
	    i__2 = k + j * b_dim1;
#line 235 "cgtsv.f"
	    i__3 = k + j * b_dim1;
#line 235 "cgtsv.f"
	    i__4 = k;
#line 235 "cgtsv.f"
	    i__5 = k + 1 + j * b_dim1;
#line 235 "cgtsv.f"
	    z__4.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, z__4.i =
		     du[i__4].r * b[i__5].i + du[i__4].i * b[i__5].r;
#line 235 "cgtsv.f"
	    z__3.r = b[i__3].r - z__4.r, z__3.i = b[i__3].i - z__4.i;
#line 235 "cgtsv.f"
	    i__6 = k;
#line 235 "cgtsv.f"
	    i__7 = k + 2 + j * b_dim1;
#line 235 "cgtsv.f"
	    z__5.r = dl[i__6].r * b[i__7].r - dl[i__6].i * b[i__7].i, z__5.i =
		     dl[i__6].r * b[i__7].i + dl[i__6].i * b[i__7].r;
#line 235 "cgtsv.f"
	    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
#line 235 "cgtsv.f"
	    z_div(&z__1, &z__2, &d__[k]);
#line 235 "cgtsv.f"
	    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 237 "cgtsv.f"
/* L40: */
#line 237 "cgtsv.f"
	}
#line 238 "cgtsv.f"
/* L50: */
#line 238 "cgtsv.f"
    }

#line 240 "cgtsv.f"
    return 0;

/*     End of CGTSV */

} /* cgtsv_ */

