#line 1 "cptts2.f"
/* cptts2.f -- translated by f2c (version 20100827).
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

#line 1 "cptts2.f"
/* > \brief \b CPTTS2 solves a tridiagonal system of the form AX=B using the L D LH factorization computed by 
spttrf. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPTTS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptts2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptts2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptts2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPTTS2( IUPLO, N, NRHS, D, E, B, LDB ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IUPLO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ) */
/*       COMPLEX            B( LDB, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPTTS2 solves a tridiagonal system of the form */
/* >    A * X = B */
/* > using the factorization A = U**H*D*U or A = L*D*L**H computed by CPTTRF. */
/* > D is a diagonal matrix specified in the vector D, U (or L) is a unit */
/* > bidiagonal matrix whose superdiagonal (subdiagonal) is specified in */
/* > the vector E, and X and B are N by NRHS matrices. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] IUPLO */
/* > \verbatim */
/* >          IUPLO is INTEGER */
/* >          Specifies the form of the factorization and whether the */
/* >          vector E is the superdiagonal of the upper bidiagonal factor */
/* >          U or the subdiagonal of the lower bidiagonal factor L. */
/* >          = 1:  A = U**H *D*U, E is the superdiagonal of U */
/* >          = 0:  A = L*D*L**H, E is the subdiagonal of L */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the tridiagonal matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the diagonal matrix D from the */
/* >          factorization A = U**H *D*U or A = L*D*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N-1) */
/* >          If IUPLO = 1, the (n-1) superdiagonal elements of the unit */
/* >          bidiagonal factor U from the factorization A = U**H*D*U. */
/* >          If IUPLO = 0, the (n-1) subdiagonal elements of the unit */
/* >          bidiagonal factor L from the factorization A = L*D*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side vectors B for the system of */
/* >          linear equations. */
/* >          On exit, the solution vectors, X. */
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

/* > \ingroup complexPTcomputational */

/*  ===================================================================== */
/* Subroutine */ int cptts2_(integer *iuplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 144 "cptts2.f"
    /* Parameter adjustments */
#line 144 "cptts2.f"
    --d__;
#line 144 "cptts2.f"
    --e;
#line 144 "cptts2.f"
    b_dim1 = *ldb;
#line 144 "cptts2.f"
    b_offset = 1 + b_dim1;
#line 144 "cptts2.f"
    b -= b_offset;
#line 144 "cptts2.f"

#line 144 "cptts2.f"
    /* Function Body */
#line 144 "cptts2.f"
    if (*n <= 1) {
#line 145 "cptts2.f"
	if (*n == 1) {
#line 145 "cptts2.f"
	    d__1 = 1. / d__[1];
#line 145 "cptts2.f"
	    csscal_(nrhs, &d__1, &b[b_offset], ldb);
#line 145 "cptts2.f"
	}
#line 147 "cptts2.f"
	return 0;
#line 148 "cptts2.f"
    }

#line 150 "cptts2.f"
    if (*iuplo == 1) {

/*        Solve A * X = B using the factorization A = U**H *D*U, */
/*        overwriting each right hand side vector with its solution. */

#line 155 "cptts2.f"
	if (*nrhs <= 2) {
#line 156 "cptts2.f"
	    j = 1;
#line 157 "cptts2.f"
L5:

/*           Solve U**H * x = b. */

#line 161 "cptts2.f"
	    i__1 = *n;
#line 161 "cptts2.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 162 "cptts2.f"
		i__2 = i__ + j * b_dim1;
#line 162 "cptts2.f"
		i__3 = i__ + j * b_dim1;
#line 162 "cptts2.f"
		i__4 = i__ - 1 + j * b_dim1;
#line 162 "cptts2.f"
		d_cnjg(&z__3, &e[i__ - 1]);
#line 162 "cptts2.f"
		z__2.r = b[i__4].r * z__3.r - b[i__4].i * z__3.i, z__2.i = b[
			i__4].r * z__3.i + b[i__4].i * z__3.r;
#line 162 "cptts2.f"
		z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
#line 162 "cptts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 163 "cptts2.f"
/* L10: */
#line 163 "cptts2.f"
	    }

/*           Solve D * U * x = b. */

#line 167 "cptts2.f"
	    i__1 = *n;
#line 167 "cptts2.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 168 "cptts2.f"
		i__2 = i__ + j * b_dim1;
#line 168 "cptts2.f"
		i__3 = i__ + j * b_dim1;
#line 168 "cptts2.f"
		i__4 = i__;
#line 168 "cptts2.f"
		z__1.r = b[i__3].r / d__[i__4], z__1.i = b[i__3].i / d__[i__4]
			;
#line 168 "cptts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 169 "cptts2.f"
/* L20: */
#line 169 "cptts2.f"
	    }
#line 170 "cptts2.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 171 "cptts2.f"
		i__1 = i__ + j * b_dim1;
#line 171 "cptts2.f"
		i__2 = i__ + j * b_dim1;
#line 171 "cptts2.f"
		i__3 = i__ + 1 + j * b_dim1;
#line 171 "cptts2.f"
		i__4 = i__;
#line 171 "cptts2.f"
		z__2.r = b[i__3].r * e[i__4].r - b[i__3].i * e[i__4].i, 
			z__2.i = b[i__3].r * e[i__4].i + b[i__3].i * e[i__4]
			.r;
#line 171 "cptts2.f"
		z__1.r = b[i__2].r - z__2.r, z__1.i = b[i__2].i - z__2.i;
#line 171 "cptts2.f"
		b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 172 "cptts2.f"
/* L30: */
#line 172 "cptts2.f"
	    }
#line 173 "cptts2.f"
	    if (j < *nrhs) {
#line 174 "cptts2.f"
		++j;
#line 175 "cptts2.f"
		goto L5;
#line 176 "cptts2.f"
	    }
#line 177 "cptts2.f"
	} else {
#line 178 "cptts2.f"
	    i__1 = *nrhs;
#line 178 "cptts2.f"
	    for (j = 1; j <= i__1; ++j) {

/*              Solve U**H * x = b. */

#line 182 "cptts2.f"
		i__2 = *n;
#line 182 "cptts2.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 183 "cptts2.f"
		    i__3 = i__ + j * b_dim1;
#line 183 "cptts2.f"
		    i__4 = i__ + j * b_dim1;
#line 183 "cptts2.f"
		    i__5 = i__ - 1 + j * b_dim1;
#line 183 "cptts2.f"
		    d_cnjg(&z__3, &e[i__ - 1]);
#line 183 "cptts2.f"
		    z__2.r = b[i__5].r * z__3.r - b[i__5].i * z__3.i, z__2.i =
			     b[i__5].r * z__3.i + b[i__5].i * z__3.r;
#line 183 "cptts2.f"
		    z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4].i - z__2.i;
#line 183 "cptts2.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 184 "cptts2.f"
/* L40: */
#line 184 "cptts2.f"
		}

/*              Solve D * U * x = b. */

#line 188 "cptts2.f"
		i__2 = *n + j * b_dim1;
#line 188 "cptts2.f"
		i__3 = *n + j * b_dim1;
#line 188 "cptts2.f"
		i__4 = *n;
#line 188 "cptts2.f"
		z__1.r = b[i__3].r / d__[i__4], z__1.i = b[i__3].i / d__[i__4]
			;
#line 188 "cptts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 189 "cptts2.f"
		for (i__ = *n - 1; i__ >= 1; --i__) {
#line 190 "cptts2.f"
		    i__2 = i__ + j * b_dim1;
#line 190 "cptts2.f"
		    i__3 = i__ + j * b_dim1;
#line 190 "cptts2.f"
		    i__4 = i__;
#line 190 "cptts2.f"
		    z__2.r = b[i__3].r / d__[i__4], z__2.i = b[i__3].i / d__[
			    i__4];
#line 190 "cptts2.f"
		    i__5 = i__ + 1 + j * b_dim1;
#line 190 "cptts2.f"
		    i__6 = i__;
#line 190 "cptts2.f"
		    z__3.r = b[i__5].r * e[i__6].r - b[i__5].i * e[i__6].i, 
			    z__3.i = b[i__5].r * e[i__6].i + b[i__5].i * e[
			    i__6].r;
#line 190 "cptts2.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 190 "cptts2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 191 "cptts2.f"
/* L50: */
#line 191 "cptts2.f"
		}
#line 192 "cptts2.f"
/* L60: */
#line 192 "cptts2.f"
	    }
#line 193 "cptts2.f"
	}
#line 194 "cptts2.f"
    } else {

/*        Solve A * X = B using the factorization A = L*D*L**H, */
/*        overwriting each right hand side vector with its solution. */

#line 199 "cptts2.f"
	if (*nrhs <= 2) {
#line 200 "cptts2.f"
	    j = 1;
#line 201 "cptts2.f"
L65:

/*           Solve L * x = b. */

#line 205 "cptts2.f"
	    i__1 = *n;
#line 205 "cptts2.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 206 "cptts2.f"
		i__2 = i__ + j * b_dim1;
#line 206 "cptts2.f"
		i__3 = i__ + j * b_dim1;
#line 206 "cptts2.f"
		i__4 = i__ - 1 + j * b_dim1;
#line 206 "cptts2.f"
		i__5 = i__ - 1;
#line 206 "cptts2.f"
		z__2.r = b[i__4].r * e[i__5].r - b[i__4].i * e[i__5].i, 
			z__2.i = b[i__4].r * e[i__5].i + b[i__4].i * e[i__5]
			.r;
#line 206 "cptts2.f"
		z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
#line 206 "cptts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 207 "cptts2.f"
/* L70: */
#line 207 "cptts2.f"
	    }

/*           Solve D * L**H * x = b. */

#line 211 "cptts2.f"
	    i__1 = *n;
#line 211 "cptts2.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 212 "cptts2.f"
		i__2 = i__ + j * b_dim1;
#line 212 "cptts2.f"
		i__3 = i__ + j * b_dim1;
#line 212 "cptts2.f"
		i__4 = i__;
#line 212 "cptts2.f"
		z__1.r = b[i__3].r / d__[i__4], z__1.i = b[i__3].i / d__[i__4]
			;
#line 212 "cptts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 213 "cptts2.f"
/* L80: */
#line 213 "cptts2.f"
	    }
#line 214 "cptts2.f"
	    for (i__ = *n - 1; i__ >= 1; --i__) {
#line 215 "cptts2.f"
		i__1 = i__ + j * b_dim1;
#line 215 "cptts2.f"
		i__2 = i__ + j * b_dim1;
#line 215 "cptts2.f"
		i__3 = i__ + 1 + j * b_dim1;
#line 215 "cptts2.f"
		d_cnjg(&z__3, &e[i__]);
#line 215 "cptts2.f"
		z__2.r = b[i__3].r * z__3.r - b[i__3].i * z__3.i, z__2.i = b[
			i__3].r * z__3.i + b[i__3].i * z__3.r;
#line 215 "cptts2.f"
		z__1.r = b[i__2].r - z__2.r, z__1.i = b[i__2].i - z__2.i;
#line 215 "cptts2.f"
		b[i__1].r = z__1.r, b[i__1].i = z__1.i;
#line 216 "cptts2.f"
/* L90: */
#line 216 "cptts2.f"
	    }
#line 217 "cptts2.f"
	    if (j < *nrhs) {
#line 218 "cptts2.f"
		++j;
#line 219 "cptts2.f"
		goto L65;
#line 220 "cptts2.f"
	    }
#line 221 "cptts2.f"
	} else {
#line 222 "cptts2.f"
	    i__1 = *nrhs;
#line 222 "cptts2.f"
	    for (j = 1; j <= i__1; ++j) {

/*              Solve L * x = b. */

#line 226 "cptts2.f"
		i__2 = *n;
#line 226 "cptts2.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 227 "cptts2.f"
		    i__3 = i__ + j * b_dim1;
#line 227 "cptts2.f"
		    i__4 = i__ + j * b_dim1;
#line 227 "cptts2.f"
		    i__5 = i__ - 1 + j * b_dim1;
#line 227 "cptts2.f"
		    i__6 = i__ - 1;
#line 227 "cptts2.f"
		    z__2.r = b[i__5].r * e[i__6].r - b[i__5].i * e[i__6].i, 
			    z__2.i = b[i__5].r * e[i__6].i + b[i__5].i * e[
			    i__6].r;
#line 227 "cptts2.f"
		    z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4].i - z__2.i;
#line 227 "cptts2.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 228 "cptts2.f"
/* L100: */
#line 228 "cptts2.f"
		}

/*              Solve D * L**H * x = b. */

#line 232 "cptts2.f"
		i__2 = *n + j * b_dim1;
#line 232 "cptts2.f"
		i__3 = *n + j * b_dim1;
#line 232 "cptts2.f"
		i__4 = *n;
#line 232 "cptts2.f"
		z__1.r = b[i__3].r / d__[i__4], z__1.i = b[i__3].i / d__[i__4]
			;
#line 232 "cptts2.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 233 "cptts2.f"
		for (i__ = *n - 1; i__ >= 1; --i__) {
#line 234 "cptts2.f"
		    i__2 = i__ + j * b_dim1;
#line 234 "cptts2.f"
		    i__3 = i__ + j * b_dim1;
#line 234 "cptts2.f"
		    i__4 = i__;
#line 234 "cptts2.f"
		    z__2.r = b[i__3].r / d__[i__4], z__2.i = b[i__3].i / d__[
			    i__4];
#line 234 "cptts2.f"
		    i__5 = i__ + 1 + j * b_dim1;
#line 234 "cptts2.f"
		    d_cnjg(&z__4, &e[i__]);
#line 234 "cptts2.f"
		    z__3.r = b[i__5].r * z__4.r - b[i__5].i * z__4.i, z__3.i =
			     b[i__5].r * z__4.i + b[i__5].i * z__4.r;
#line 234 "cptts2.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 234 "cptts2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 236 "cptts2.f"
/* L110: */
#line 236 "cptts2.f"
		}
#line 237 "cptts2.f"
/* L120: */
#line 237 "cptts2.f"
	    }
#line 238 "cptts2.f"
	}
#line 239 "cptts2.f"
    }

#line 241 "cptts2.f"
    return 0;

/*     End of CPTTS2 */

} /* cptts2_ */

