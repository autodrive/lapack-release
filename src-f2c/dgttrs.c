#line 1 "dgttrs.f"
/* dgttrs.f -- translated by f2c (version 20100827).
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

#line 1 "dgttrs.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b DGTTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGTTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgttrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgttrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgttrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGTTRS solves one of the systems of equations */
/* >    A*X = B  or  A**T*X = B, */
/* > with a tridiagonal matrix A using the LU factorization computed */
/* > by DGTTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations. */
/* >          = 'N':  A * X = B  (No transpose) */
/* >          = 'T':  A**T* X = B  (Transpose) */
/* >          = 'C':  A**T* X = B  (Conjugate transpose = Transpose) */
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
/* >          DL is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) elements of the first super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is DOUBLE PRECISION array, dimension (N-2) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the matrix of right hand side vectors B. */
/* >          On exit, B is overwritten by the solution vectors X. */
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
/* Subroutine */ int dgttrs_(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, 
	integer *ipiv, doublereal *b, integer *ldb, integer *info, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, jb, nb;
    extern /* Subroutine */ int dgtts2_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer itrans;
    static logical notran;


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

#line 173 "dgttrs.f"
    /* Parameter adjustments */
#line 173 "dgttrs.f"
    --dl;
#line 173 "dgttrs.f"
    --d__;
#line 173 "dgttrs.f"
    --du;
#line 173 "dgttrs.f"
    --du2;
#line 173 "dgttrs.f"
    --ipiv;
#line 173 "dgttrs.f"
    b_dim1 = *ldb;
#line 173 "dgttrs.f"
    b_offset = 1 + b_dim1;
#line 173 "dgttrs.f"
    b -= b_offset;
#line 173 "dgttrs.f"

#line 173 "dgttrs.f"
    /* Function Body */
#line 173 "dgttrs.f"
    *info = 0;
#line 174 "dgttrs.f"
    notran = *(unsigned char *)trans == 'N' || *(unsigned char *)trans == 'n';
#line 175 "dgttrs.f"
    if (! notran && ! (*(unsigned char *)trans == 'T' || *(unsigned char *)
	    trans == 't') && ! (*(unsigned char *)trans == 'C' || *(unsigned 
	    char *)trans == 'c')) {
#line 177 "dgttrs.f"
	*info = -1;
#line 178 "dgttrs.f"
    } else if (*n < 0) {
#line 179 "dgttrs.f"
	*info = -2;
#line 180 "dgttrs.f"
    } else if (*nrhs < 0) {
#line 181 "dgttrs.f"
	*info = -3;
#line 182 "dgttrs.f"
    } else if (*ldb < max(*n,1)) {
#line 183 "dgttrs.f"
	*info = -10;
#line 184 "dgttrs.f"
    }
#line 185 "dgttrs.f"
    if (*info != 0) {
#line 186 "dgttrs.f"
	i__1 = -(*info);
#line 186 "dgttrs.f"
	xerbla_("DGTTRS", &i__1, (ftnlen)6);
#line 187 "dgttrs.f"
	return 0;
#line 188 "dgttrs.f"
    }

/*     Quick return if possible */

#line 192 "dgttrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 192 "dgttrs.f"
	return 0;
#line 192 "dgttrs.f"
    }

/*     Decode TRANS */

#line 197 "dgttrs.f"
    if (notran) {
#line 198 "dgttrs.f"
	itrans = 0;
#line 199 "dgttrs.f"
    } else {
#line 200 "dgttrs.f"
	itrans = 1;
#line 201 "dgttrs.f"
    }

/*     Determine the number of right-hand sides to solve at a time. */

#line 205 "dgttrs.f"
    if (*nrhs == 1) {
#line 206 "dgttrs.f"
	nb = 1;
#line 207 "dgttrs.f"
    } else {
/* Computing MAX */
#line 208 "dgttrs.f"
	i__1 = 1, i__2 = ilaenv_(&c__1, "DGTTRS", trans, n, nrhs, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
#line 208 "dgttrs.f"
	nb = max(i__1,i__2);
#line 209 "dgttrs.f"
    }

#line 211 "dgttrs.f"
    if (nb >= *nrhs) {
#line 212 "dgttrs.f"
	dgtts2_(&itrans, n, nrhs, &dl[1], &d__[1], &du[1], &du2[1], &ipiv[1], 
		&b[b_offset], ldb);
#line 213 "dgttrs.f"
    } else {
#line 214 "dgttrs.f"
	i__1 = *nrhs;
#line 214 "dgttrs.f"
	i__2 = nb;
#line 214 "dgttrs.f"
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
#line 215 "dgttrs.f"
	    i__3 = *nrhs - j + 1;
#line 215 "dgttrs.f"
	    jb = min(i__3,nb);
#line 216 "dgttrs.f"
	    dgtts2_(&itrans, n, &jb, &dl[1], &d__[1], &du[1], &du2[1], &ipiv[
		    1], &b[j * b_dim1 + 1], ldb);
#line 218 "dgttrs.f"
/* L10: */
#line 218 "dgttrs.f"
	}
#line 219 "dgttrs.f"
    }

/*     End of DGTTRS */

#line 223 "dgttrs.f"
    return 0;
} /* dgttrs_ */

