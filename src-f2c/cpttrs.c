#line 1 "cpttrs.f"
/* cpttrs.f -- translated by f2c (version 20100827).
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

#line 1 "cpttrs.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CPTTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPTTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpttrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpttrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpttrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS */
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
/* > CPTTRS solves a tridiagonal system of the form */
/* >    A * X = B */
/* > using the factorization A = U**H*D*U or A = L*D*L**H computed by CPTTRF. */
/* > D is a diagonal matrix specified in the vector D, U (or L) is a unit */
/* > bidiagonal matrix whose superdiagonal (subdiagonal) is specified in */
/* > the vector E, and X and B are N by NRHS matrices. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies the form of the factorization and whether the */
/* >          vector E is the superdiagonal of the upper bidiagonal factor */
/* >          U or the subdiagonal of the lower bidiagonal factor L. */
/* >          = 'U':  A = U**H*D*U, E is the superdiagonal of U */
/* >          = 'L':  A = L*D*L**H, E is the subdiagonal of L */
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
/* >          factorization A = U**H*D*U or A = L*D*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N-1) */
/* >          If UPLO = 'U', the (n-1) superdiagonal elements of the unit */
/* >          bidiagonal factor U from the factorization A = U**H*D*U. */
/* >          If UPLO = 'L', the (n-1) subdiagonal elements of the unit */
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
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
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
/* Subroutine */ int cpttrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, jb, nb, iuplo;
    static logical upper;
    extern /* Subroutine */ int cptts2_(integer *, integer *, integer *, 
	    doublereal *, doublecomplex *, doublecomplex *, integer *), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 158 "cpttrs.f"
    /* Parameter adjustments */
#line 158 "cpttrs.f"
    --d__;
#line 158 "cpttrs.f"
    --e;
#line 158 "cpttrs.f"
    b_dim1 = *ldb;
#line 158 "cpttrs.f"
    b_offset = 1 + b_dim1;
#line 158 "cpttrs.f"
    b -= b_offset;
#line 158 "cpttrs.f"

#line 158 "cpttrs.f"
    /* Function Body */
#line 158 "cpttrs.f"
    *info = 0;
#line 159 "cpttrs.f"
    upper = *(unsigned char *)uplo == 'U' || *(unsigned char *)uplo == 'u';
#line 160 "cpttrs.f"
    if (! upper && ! (*(unsigned char *)uplo == 'L' || *(unsigned char *)uplo 
	    == 'l')) {
#line 161 "cpttrs.f"
	*info = -1;
#line 162 "cpttrs.f"
    } else if (*n < 0) {
#line 163 "cpttrs.f"
	*info = -2;
#line 164 "cpttrs.f"
    } else if (*nrhs < 0) {
#line 165 "cpttrs.f"
	*info = -3;
#line 166 "cpttrs.f"
    } else if (*ldb < max(1,*n)) {
#line 167 "cpttrs.f"
	*info = -7;
#line 168 "cpttrs.f"
    }
#line 169 "cpttrs.f"
    if (*info != 0) {
#line 170 "cpttrs.f"
	i__1 = -(*info);
#line 170 "cpttrs.f"
	xerbla_("CPTTRS", &i__1, (ftnlen)6);
#line 171 "cpttrs.f"
	return 0;
#line 172 "cpttrs.f"
    }

/*     Quick return if possible */

#line 176 "cpttrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 176 "cpttrs.f"
	return 0;
#line 176 "cpttrs.f"
    }

/*     Determine the number of right-hand sides to solve at a time. */

#line 181 "cpttrs.f"
    if (*nrhs == 1) {
#line 182 "cpttrs.f"
	nb = 1;
#line 183 "cpttrs.f"
    } else {
/* Computing MAX */
#line 184 "cpttrs.f"
	i__1 = 1, i__2 = ilaenv_(&c__1, "CPTTRS", uplo, n, nrhs, &c_n1, &c_n1,
		 (ftnlen)6, (ftnlen)1);
#line 184 "cpttrs.f"
	nb = max(i__1,i__2);
#line 185 "cpttrs.f"
    }

/*     Decode UPLO */

#line 189 "cpttrs.f"
    if (upper) {
#line 190 "cpttrs.f"
	iuplo = 1;
#line 191 "cpttrs.f"
    } else {
#line 192 "cpttrs.f"
	iuplo = 0;
#line 193 "cpttrs.f"
    }

#line 195 "cpttrs.f"
    if (nb >= *nrhs) {
#line 196 "cpttrs.f"
	cptts2_(&iuplo, n, nrhs, &d__[1], &e[1], &b[b_offset], ldb);
#line 197 "cpttrs.f"
    } else {
#line 198 "cpttrs.f"
	i__1 = *nrhs;
#line 198 "cpttrs.f"
	i__2 = nb;
#line 198 "cpttrs.f"
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 199 "cpttrs.f"
	    i__3 = *nrhs - j + 1;
#line 199 "cpttrs.f"
	    jb = min(i__3,nb);
#line 200 "cpttrs.f"
	    cptts2_(&iuplo, n, &jb, &d__[1], &e[1], &b[j * b_dim1 + 1], ldb);
#line 201 "cpttrs.f"
/* L10: */
#line 201 "cpttrs.f"
	}
#line 202 "cpttrs.f"
    }

#line 204 "cpttrs.f"
    return 0;

/*     End of CPTTRS */

} /* cpttrs_ */

